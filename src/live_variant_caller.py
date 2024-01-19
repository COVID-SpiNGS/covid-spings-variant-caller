import functools
import operator
import math
from typing import List

import pysam
import numpy as np

from tqdm import tqdm
from pysam import AlignedSegment
from time import strftime, localtime

from variant_caller.models import Site, Variant
import variant_caller.utils as u
import variant_caller.vcf_file_constants as c
import file_util as fu

class LiveVariantCaller:
    def __init__(self, reference_fasta: str, min_base_quality: int, min_mapping_quality: int, min_total_depth: int,
                 min_allele_depth: int, min_evidence_ratio: float, max_variants: int):
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality
        self.min_total_depth = min_total_depth
        self.min_allele_depth = min_allele_depth
        self.min_evidence_ratio = min_evidence_ratio
        self.max_variants = max_variants
        self.fasta_file = pysam.FastaFile(reference_fasta)
        self.memory = {}
        self.reset_memory()

    def __del__(self):
        self.fasta_file.close()

    def reset_memory(self):
        self.memory: dict[int, Site] = {}


    def process_bam(self, input_bam: str, reference_index=0):
        bam_file = pysam.AlignmentFile(input_bam, 'rb')
        pileup_columns =  bam_file.pileup(
            min_mapping_quality=self.min_mapping_quality,
            min_base_quality=self.min_base_quality,
            reference=self.fasta_file.references[reference_index]
        )

        timestamp = strftime('[%Y-%m-%d %H:%M:%S]', localtime())
        progress_bar = tqdm(
            pileup_columns,
            desc=f'{timestamp} Processing {input_bam}',
            total= bam_file.get_reference_length(self.fasta_file.references[reference_index])
        )

        for pileup_column in progress_bar:
            self.process_pileup_column(pileup_column)

        bam_file.close()

    def process_pileup_column(self, pileup_column: AlignedSegment):
        total_depth = len(pileup_column.pileups)

        if pileup_column.reference_pos not in self.memory:
            reference = self.fasta_file.fetch(reference=pileup_column.reference_name)

            self.memory[pileup_column.reference_pos] = {
                c.VCF_REFERENCE: reference[pileup_column.reference_pos],
                c.VCF_TOTAL_DEPTH_KEY: total_depth,
                c.VCF_SNVS: {},
                c.VCF_INDELS: {}
            }
        else:
            self.memory[pileup_column.reference_pos][c.VCF_TOTAL_DEPTH_KEY] += total_depth

        for pileup in pileup_column.pileups:
            self.process_pileup_at_position(pileup_column.reference_pos, pileup)

    def process_pileup_at_position(self, position: int, pileup):
        self.process_snv(position, pileup)
        # self.process_indel(position, pileup)

    def process_snv(self, position, pileup):
        if not pileup.is_del and not pileup.is_refskip:
            snv = pileup.alignment.query_sequence[pileup.query_position]

            if snv not in self.memory[position][c.VCF_SNVS].keys():
                self.memory[position][c.VCF_SNVS][snv] = []

            self.memory[position][c.VCF_SNVS][snv].append(pileup.alignment.query_qualities[pileup.query_position])

    def process_indel(self, position, pileup):
        if pileup.is_del or pileup.is_refskip:
            indel = '-' if pileup.is_del else f'+{pileup.alignment.query_sequence[pileup.query_position]}'

            if indel not in self.memory[position][c.VCF_INDELS].keys():
                # We could store more information here in the memory. But as the Base Qualty is the the only information that matters for further calculation we 
                # save memory and only store them

                self.memory[position][c.VCF_INDELS][indel] = []

            if pileup.is_refskip:
                self.memory[position][c.VCF_INDELS][indel].append(pileup.alignment.query_qualities[pileup.query_position])
            else:
                self.memory[position][c.VCF_INDELS][indel].append(None)

    def prepare_variants(self):
        timestamp = strftime('[%Y-%m-%d %H:%M:%S]', localtime())
        progress_bar = tqdm(
            self.memory,
            desc=f'{timestamp} Calculating statistics',
            total=len(self.memory.keys())
        )

        variants: List[Variant] = []

        for position in progress_bar:
            if self.memory[position][c.VCF_TOTAL_DEPTH_KEY] >= self.min_total_depth:
                snvs = {
                    allele: [
                        u.from_phred_score(quality)[1]
                        for quality in self.memory[position][c.VCF_SNVS][allele]
                    ]
                    for allele in self.memory[position][c.VCF_SNVS].keys()
                }

                snvs2 = {
                    allele: [
                        u.from_phred_score(quality)[0]
                        for quality in self.memory[position][c.VCF_SNVS][allele]
                    ]
                    for allele in self.memory[position][c.VCF_SNVS].keys()
                }

                snvs_tuples = [(key, value) for key in snvs2 for value in snvs2[key]]
            
                # MAGIC HAPPENS HERE 

                if len(snvs_tuples) > 0:
                    x = u.get_likelihood(snvs_tuples)

                    genotype_likelihoods = u.extract_base_likelihood(x, snvs_tuples, snvs)
                
                    print(genotype_likelihoods)
                #genotype_likelihoodsss = {
                #    allele: u.genotype_likelihood2(allele, snvs)
                #   for allele in snvs.keys()
                #}

               
                    sum_genotype_likelihoods = functools.reduce(operator.add, genotype_likelihoods.values(), 0.0)
                    sum_genotype_likelihoods = sum_genotype_likelihoods if sum_genotype_likelihoods != 0 else 1.0

                for allele in snvs.keys():
                    allele_depth = len(snvs[allele])

                    filter_constrains = [
                        self.memory[position][c.VCF_REFERENCE] != allele,
                        allele_depth >= self.min_allele_depth,
                        allele_depth / self.memory[position][c.VCF_TOTAL_DEPTH_KEY] >= self.min_evidence_ratio
                    ]

                    if all(filter_constrains):
                        genotype_likelihood = genotype_likelihoods[allele]

                        if genotype_likelihood != 0:
                            gl = math.log10(genotype_likelihood)
                            pl = round(-10.0 * gl)
                        else:
                            gl = 0
                            pl = 0

                        # MAGIC HAPPENS HERE TOO I GUESS ?
                        score = u.to_phred_score(1.0 - (genotype_likelihoods[allele] / sum_genotype_likelihoods))
                        qual = np.mean(snvs[allele])

                        variants.append({
                            c.VCF_START: position,
                            c.VCF_STOP: position + 1,
                            c.VCF_ALLELES: (
                                self.memory[position][c.VCF_REFERENCE],
                                allele
                            ),
                            c.VCF_QUAL: qual,
                            c.VCF_INFO: {
                                c.VCF_DP: self.memory[position][c.VCF_TOTAL_DEPTH_KEY],
                                c.VCF_AD: allele_depth,
                                c.VCF_GL: gl,
                                c.VCF_PL: pl,
                                c.VCF_SCORE: score
                            }
                        })

                for indel in self.memory[position][c.VCF_INDELS].keys():
                    allele_depth = len(self.memory[position][c.VCF_INDELS][indel])

                    filter_constrains = [
                        allele_depth >= self.min_allele_depth,
                        allele_depth / self.memory[position][c.VCF_TOTAL_DEPTH_KEY] >= self.min_evidence_ratio
                    ]

                    if all(filter_constrains):
                        if indel == '-':
                            variants.append({
                                c.VCF_START: position,
                                c.VCF_STOP: position + 1,
                                c.VCF_ALLELES: (
                                    self.memory[position][c.VCF_REFERENCE],
                                    '*'
                                ),
                                c.VCF_QUAL: 0,
                                c.VCF_INFO: {
                                    c.VCF_DP: self.memory[position][c.VCF_TOTAL_DEPTH_KEY],
                                    c.VCF_AD: allele_depth,
                                    c.VCF_GL: 0,
                                    c.VCF_PL: 0,
                                    c.VCF_SCORE: 0
                                }
                            })
                        else:
                            variants.append({
                                c.VCF_START: position,
                                c.VCF_STOP: position + 1,
                                c.VCF_ALLELES: (
                                    '*',
                                    indel[1:]
                                ),
                                c.VCF_QUAL: 0,
                                c.VCF_INFO: {
                                    c.VCF_DP: self.memory[position][c.VCF_TOTAL_DEPTH_KEY],
                                    c.VCF_ED: allele_depth,
                                    c.VCF_GL: 0,
                                    c.VCF_PL: 0,
                                    c.VCF_SCORE: 0
                                }
                            })

        return variants


    def prev_variant(self, variants: List[Variant], variant: Variant):
        return next(
            (
                v for v in variants
                if v[c.VCF_START] == variant[c.VCF_START] - 1
            ),
            None
        )

    def next_variant(self, variants: List[Variant], variant: Variant):
        return next(
            (
                v for v in variants
                if v[c.VCF_START] == variant[c.VCF_START] + 1
            ),
            None
        )

    def concat_deletions(self, variants: List[Variant]):
        concatinated_variants: List[Variant] = []
        current_variant: Variant = None

        for variant in variants:
            if variant[c.VCF_ALLELES][1] == '*':
                next_variant = self.next_variant(variants, variant)

                if next_variant:
                    if not current_variant:
                        current_variant = variant
                    else:
                        current_variant = {
                            c.VCF_START: current_variant[c.VCF_START],
                            c.VCF_STOP: variant[c.VCF_STOP],
                            c.VCF_ALLELES: (
                                f'{current_variant[c.VCF_ALLELES][0]}{variant[c.VCF_ALLELES][0]}',
                                '*'
                            ),
                            c.VCF_QUAL: variant[c.VCF_QUAL],  # must be combined
                            c.VCF_INFO: variant[c.VCF_INFO]  # must be combined

                        }
                else:
                    if current_variant:
                        concatinated_variants.append(current_variant)
                        current_variant = None


            else:
                concatinated_variants.append(variant)

        return concatinated_variants

    def concat_insertions(self, variants: List[Variant]):
        return variants
    
    def write_vcf(self, output_vcf: str):
            print("VCF output", output_vcf)
            vcf_header = pysam.VariantHeader()
    
            vcf_header.add_meta(c.VCF_INFO, items=[
                (c.VCF_ID, c.VCF_DP),
                (c.VCF_NUMBER, 1),
                (c.VCF_TYPE, c.VCF_TYPE_INTEGER),
                (c.VCF_DESCRIPTION, c.VCF_TOTAL_DEPTH_STR)
            ])
    
            vcf_header.add_meta(c.VCF_INFO, items=[
                (c.VCF_ID, c.VCF_AD),
                (c.VCF_NUMBER, 1),
                (c.VCF_TYPE, c.VCF_TYPE_INTEGER),
                (c.VCF_DESCRIPTION, c.VCF_ALLELE_DEPTH)
            ])
    
            vcf_header.add_meta(c.VCF_INFO, items=[
                (c.VCF_ID, c.VCF_GL),
                (c.VCF_NUMBER, 1),
                (c.VCF_TYPE, c.VCF_TYPE_FLOAT),
                (c.VCF_DESCRIPTION,
                 'Genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields')
            ])
    
            vcf_header.add_meta(c.VCF_INFO, items=[
                (c.VCF_ID, c.VCF_PL),
                (c.VCF_NUMBER, 1),
                (c.VCF_TYPE, c.VCF_TYPE_INTEGER),
                (c.VCF_DESCRIPTION,
                 'The phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field)')
            ])
    
            vcf_header.add_meta(c.VCF_INFO, items=[
                (c.VCF_ID, c.VCF_SCORE),
                (c.VCF_NUMBER, 1),
                (c.VCF_TYPE, c.VCF_TYPE_FLOAT),
                (c.VCF_DESCRIPTION, 'Custom scoring function')
            ])
    
            for reference in self.fasta_file.references:
                vcf_header.contigs.add(
                    reference,
                    self.fasta_file.get_reference_length(reference)
                )
    
            vcf_file = pysam.VariantFile(output_vcf, mode='w', header=vcf_header)
    
            variants = self.prepare_variants()
            # gvariants = self.concat_deletions(variants)
    
            

            for index, variant in enumerate(
                    sorted(variants, key=lambda variant: (variant[c.VCF_START], variant[c.VCF_INFO][c.VCF_SCORE]))):
                vcf_file.write(
                    vcf_file.new_record(
                        start=variant[c.VCF_START],
                        stop=variant[c.VCF_STOP],
                        alleles=variant[c.VCF_ALLELES],
                        qual=variant[c.VCF_QUAL],
                        info=variant[c.VCF_INFO],
                    )
                )
    
            vcf_file.close()
    