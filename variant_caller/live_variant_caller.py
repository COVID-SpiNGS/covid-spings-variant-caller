import functools
import operator
import pickle
import math
from typing import List

import os

import pysam
import numpy as np

from tqdm import tqdm
from pysam import AlignedSegment
from time import strftime, localtime

from .structs import Site, Variant
from .utils import genotype_likelihood, from_phred_scale, to_phred_scale



class LiveVariantCaller:
    def __init__(self, referenceFasta: str, minBaseQuality: int, minMappingQuality: int, minTotalDepth: int,
                 minAlleleDepth: int, minEvidenceRatio: float, maxVariants: int):
        self.minBaseQuality = minBaseQuality
        self.minMappingQuality = minMappingQuality
        self.minTotalDepth = minTotalDepth
        self.minAlleleDepth = minAlleleDepth
        self.minEvidenceRatio = minEvidenceRatio
        self.maxVariants = maxVariants
        self.fastaFile = pysam.FastaFile(referenceFasta)
        self.memory = {}
        self.reset_memory()

    def __del__(self):
        self.fastaFile.close()

    def reset_memory(self):
        self.memory: dict[int, Site] = {}

    def create_checkpoint(self, filename):
        timestamp = strftime('[%Y-%m-%d %H:%M:%S]', localtime())
        print(f'{timestamp} Creating checkpoint {filename}')

        file = open(filename, 'wb')
        pickle.dump(self.memory, file)
        file.close()

    def load_checkpoint(self, filename):
        timestamp = strftime('[%Y-%m-%d %H:%M:%S]', localtime())
        print(f'{timestamp} Loading checkpoint {filename}')

        file = open(filename, 'rb')
        self.memory = pickle.load(file)
        file.close()

    def process_bam(self, inputBam: str, referenceIndex=0):
        bamFile = pysam.AlignmentFile(inputBam, 'rb')

        pileupColumns = bamFile.pileup(
            min_mapping_quality=self.minMappingQuality,
            min_base_quality=self.minBaseQuality,
            reference=self.fastaFile.references[referenceIndex]
        )

        timestamp = strftime('[%Y-%m-%d %H:%M:%S]', localtime())
        progressBar = tqdm(
            pileupColumns,
            desc=f'{timestamp} Processing {inputBam}',
            total=bamFile.get_reference_length(self.fastaFile.references[referenceIndex])
        )

        for pileupColumn in progressBar:
            self.process_pileup_column(pileupColumn)

        bamFile.close()

    def process_pileup_column(self, pileupColumn: AlignedSegment):
        totalDepth = len(pileupColumn.pileups)

        if pileupColumn.reference_pos not in self.memory:
            reference = self.fastaFile.fetch(reference=pileupColumn.reference_name)

            self.memory[pileupColumn.reference_pos] = {
                'reference': reference[pileupColumn.reference_pos],
                'totalDepth': totalDepth,
                'snvs': {},
                'indels': {}
            }
        else:
            self.memory[pileupColumn.reference_pos]['totalDepth'] += totalDepth

        for pileup in pileupColumn.pileups:
            self.process_pileup_at_position(pileupColumn.reference_pos, pileup)

    def process_pileup_at_position(self, position: int, pileup):
        self.process_svn(position, pileup)
        # self.process_indel(position, pileup)

    def process_svn(self, position, pileup):
        if not pileup.is_del and not pileup.is_refskip:
            svn = pileup.alignment.query_sequence[pileup.query_position]

            if svn not in self.memory[position]['snvs'].keys():
                self.memory[position]['snvs'][svn] = []

            self.memory[position]['snvs'][svn].append(pileup.alignment.query_qualities[pileup.query_position])

    def process_indel(self, position, pileup):
        if pileup.is_del or pileup.is_refskip:
            indel = '-' if pileup.is_del else f'+{pileup.alignment.query_sequence[pileup.query_position]}'

            if indel not in self.memory[position]['indels'].keys():
                # We could store more information here in the memory. But as the Base Qualty is the the only information that matters for further calculation we 
                # save memory and only store them

                self.memory[position]['indels'][indel] = []

            if pileup.is_refskip:
                self.memory[position]['indels'][indel].append(pileup.alignment.query_qualities[pileup.query_position])
            else:
                self.memory[position]['indels'][indel].append(None)

    def prepare_variants(self):
        timestamp = strftime('[%Y-%m-%d %H:%M:%S]', localtime())
        progressBar = tqdm(
            self.memory,
            desc=f'{timestamp} Calculating statistics',
            total=len(self.memory.keys())
        )

        variants: List[Variant] = []

        for position in progressBar:
            if self.memory[position]['totalDepth'] >= self.minTotalDepth:
                snvs = {
                    allele: [
                        from_phred_scale(quality)
                        for quality in self.memory[position]['snvs'][allele]
                    ]
                    for allele in self.memory[position]['snvs'].keys()
                }

                genotypeLikelihoods = {
                    allele: genotype_likelihood(allele, snvs)
                    for allele in snvs.keys()
                }

                sumGenotypeLikelihoods = functools.reduce(operator.add, genotypeLikelihoods.values(), 0.0)
                sumGenotypeLikelihoods = sumGenotypeLikelihoods if sumGenotypeLikelihoods != 0 else 1.0

                for allele in snvs.keys():
                    alleleDepth = len(snvs[allele])

                    filterConstrains = [
                        self.memory[position]['reference'] != allele,
                        alleleDepth >= self.minAlleleDepth,
                        alleleDepth / self.memory[position]['totalDepth'] >= self.minEvidenceRatio
                    ]

                    if all(filterConstrains):
                        genotypeLikelihood = genotypeLikelihoods[allele]

                        if genotypeLikelihood != 0:
                            gl = math.log10(genotypeLikelihood)
                            pl = round(-10.0 * gl)
                        else:
                            gl = 0
                            pl = 0

                        score = to_phred_scale(1.0 - (genotypeLikelihoods[allele] / sumGenotypeLikelihoods))
                        qual = np.mean(snvs[allele])

                        variants.append({
                            'start': position,
                            'stop': position + 1,
                            'alleles': (
                                self.memory[position]['reference'],
                                allele
                            ),
                            'qual': qual,
                            'info': {
                                'DP': self.memory[position]['totalDepth'],
                                'AD': alleleDepth,
                                'GL': gl,
                                'PL': pl,
                                'SCORE': score
                            }
                        })

                for indel in self.memory[position]['indels'].keys():
                    alleleDepth = len(self.memory[position]['indels'][indel])

                    filterConstrains = [
                        alleleDepth >= self.minAlleleDepth,
                        alleleDepth / self.memory[position]['totalDepth'] >= self.minEvidenceRatio
                    ]

                    if all(filterConstrains):
                        if indel == '-':
                            variants.append({
                                'start': position,
                                'stop': position + 1,
                                'alleles': (
                                    self.memory[position]['reference'],
                                    '*'
                                ),
                                'qual': 0,
                                'info': {
                                    'DP': self.memory[position]['totalDepth'],
                                    'AD': alleleDepth,
                                    'GL': 0,
                                    'PL': 0,
                                    'SCORE': 0
                                }
                            })
                        else:
                            variants.append({
                                'start': position,
                                'stop': position + 1,
                                'alleles': (
                                    '*',
                                    indel[1:]
                                ),
                                'qual': 0,
                                'info': {
                                    'DP': self.memory[position]['totalDepth'],
                                    'ED': alleleDepth,
                                    'GL': 0,
                                    'PL': 0,
                                    'SCORE': 0
                                }
                            })

        return variants

    def write_vcf(self, outputVfc: str):
        vcfHeader = pysam.VariantHeader()

        vcfHeader.add_meta('INFO', items=[
            ('ID', 'DP'),
            ('Number', 1),
            ('Type', 'Integer'),
            ('Description', 'Total Depth')
        ])

        vcfHeader.add_meta('INFO', items=[
            ('ID', 'AD'),
            ('Number', 1),
            ('Type', 'Integer'),
            ('Description', 'Allele Depth')
        ])

        vcfHeader.add_meta('INFO', items=[
            ('ID', 'GL'),
            ('Number', 1),
            ('Type', 'Float'),
            ('Description',
             'Genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields')
        ])

        vcfHeader.add_meta('INFO', items=[
            ('ID', 'PL'),
            ('Number', 1),
            ('Type', 'Integer'),
            ('Description',
             'The phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field)')
        ])

        vcfHeader.add_meta('INFO', items=[
            ('ID', 'SCORE'),
            ('Number', 1),
            ('Type', 'Float'),
            ('Description', 'Custom scoring function')
        ])

        for reference in self.fastaFile.references:
            vcfHeader.contigs.add(
                reference,
                self.fastaFile.get_reference_length(reference)
            )

        vcfFile = pysam.VariantFile(outputVfc, mode='w', header=vcfHeader)

        variants = self.prepare_variants()
        # gvariants = self.concat_deletions(variants)

        for index, variant in enumerate(
                sorted(variants, key=lambda variant: (variant['start'], variant['info']['SCORE']))):
            vcfFile.write(
                vcfFile.new_record(
                    start=variant['start'],
                    stop=variant['stop'],
                    alleles=variant['alleles'],
                    qual=variant['qual'],
                    info=variant['info'],
                )
            )

        vcfFile.close()

    def prev_variant(self, variants: List[Variant], variant: Variant):
        return next(
            (
                v for v in variants
                if v['start'] == variant['start'] - 1
            ),
            None
        )

    def next_variant(self, variants: List[Variant], variant: Variant):
        return next(
            (
                v for v in variants
                if v['start'] == variant['start'] + 1
            ),
            None
        )

    def concat_deletions(self, variants: List[Variant]):
        concatinatedVariants: List[Variant] = []
        currentVariant: Variant = None

        for variant in variants:
            if variant['alleles'][1] == '*':
                nextVariant = self.next_variant(variants, variant)

                if nextVariant:
                    if not currentVariant:
                        currentVariant = variant
                    else:
                        currentVariant = {
                            'start': currentVariant['start'],
                            'stop': variant['stop'],
                            'alleles': (
                                f'{currentVariant["alleles"][0]}{variant["alleles"][0]}',
                                '*'
                            ),
                            'qual': variant['qual'],  # must be combined
                            'info': variant['info']  # must be combined

                        }
                else:
                    if currentVariant:
                        concatinatedVariants.append(currentVariant)
                        currentVariant = None


            else:
                concatinatedVariants.append(variant)

        return concatinatedVariants

    def concat_insertions(self, variants: List[Variant]):
        return variants
