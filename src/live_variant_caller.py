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
import pickle
from config_util import logging as log


class LiveVariantCaller:
    """
    A class for calling genetic variants in real-time from sequencing data.

    This class provides functionalities to identify variants based on various criteria such as base quality, mapping quality, allele depth, and evidence ratio. It leverages a reference FASTA file for comparison and maintains an in-memory record of variant sites.

    Attributes:
        min_base_quality (int): Minimum base quality required for a base to be considered for variant calling.
        min_mapping_quality (int): Minimum mapping quality required for a read to be considered in variant calling.
        min_total_depth (int): Minimum total depth (number of reads) at a site to consider it for variant calling.
        min_allele_depth (int): Minimum depth of a specific allele required to consider it a variant.
        min_evidence_ratio (float): Minimum ratio of allele depth to total depth to consider a variant.
        max_variants (int): Maximum number of variants to store in memory.
        fasta_file (pysam.FastaFile): A pysam FastaFile object for accessing the reference genome.
        memory (dict[int, Site]): A dictionary to store information about variant sites, keyed by their genomic position.

    """

    def __init__(self, reference_fasta: str, min_base_quality: int, min_mapping_quality: int, min_total_depth: int,
                 min_allele_depth: int, min_evidence_ratio: float, max_variants: int):
        
        """
        Constructor that sets up the LiveVariantCaller with the necessary parameters and resources for variant calling.

        Args:
            reference_fasta (str): The file path to the reference genome FASTA file. This file is used as a baseline to compare against sequencing data for variant calling.
            min_base_quality (int): The minimum quality score required for a base to be considered reliable for variant calling. Quality scores lower than this threshold are ignored.
            min_mapping_quality (int): The minimum mapping quality score required for a read. Reads with mapping quality below this threshold are disregarded in variant calling.
            min_total_depth (int): The minimum number of reads (depth) that must cover a site for the site to be considered in variant calling.
            min_allele_depth (int): The minimum number of reads supporting an allele to consider it a variant.
            min_evidence_ratio (float): The minimum ratio of the number of reads supporting an allele to the total number of reads covering a site, required to consider it a variant.
            max_variants (int): The maximum number of variant sites to keep in memory.

        """
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
        """
        Destructor for the LiveVariantCaller class.

        The method closes the FASTA file associated with the instance, ensuring that there are no file handle leaks and that the file is properly closed when the LiveVariantCaller object is no longer in use.

        """
        self.fasta_file.close()


    
    def reset_memory(self):

        """
        Resets the memory of the LiveVariantCaller instance.

        This method clears the current memory storage used for tracking variant sites. It initializes the 'memory' attribute to an empty dictionary. This is useful for starting fresh variant calling without creating a new instance of the LiveVariantCaller.

        The 'memory' attribute holds information about variant sites. Each entry in this dictionary represents a variant site, with the key being the genomic position (integer) and the value being a 'Site' object containing variant information.

        """
        self.memory: dict[int, Site] = {}


    def process_bam(self, input_bam: str, reference_index=0):
        """
        Processes a BAM file to identify variant sites.

        This method reads a BAM file, iterates over its pileup columns, and processes each column to identify potential variant sites. The method utilizes the configured thresholds for base and mapping quality to filter reads.

        Args:
            input_bam (str): The file path to the input BAM file. This file contains the aligned sequencing reads.
            reference_index (int, optional): The index of the reference sequence in the reference FASTA file to be used for pileup. Defaults to 0, indicating the first reference sequence.

        The method uses pysam to open and iterate through the BAM file. For each pileup column, it calls the 'process_pileup_column' method to evaluate potential variants. The progress of processing is displayed using a progress bar (tqdm).

        A timestamp is added to the progress bar description to track when the processing started.
        """

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

        """
        Processes a single pileup column from a BAM file to identify and record variant information.

        This method examines each pileup column in a BAM file to evaluate the presence of potential variants (SNVs and INDELs) at that specific position in the reference sequence.

        Args:
            pileup_column (pysam.AlignedSegment): A pileup column object representing a single position in the reference sequence, along with all the reads aligned to that position.

        The method first calculates the total depth (number of reads) at the pileup column's position. It then checks if this position is already recorded in the 'memory'. If not, it initializes a dictionary at this position in 'memory' with the reference base, total depth, and placeholders for SNVs and INDELs.

        If the position is already in memory, the method updates the total depth.

        For each read in the pileup column, the method calls 'process_pileup_at_position' to evaluate individual read alignment for potential variants.

        """

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
        """
        Processes a single pileup at a given position to identify single nucleotide variants (SNVs).

        This method is a key part of the variant calling process, as it analyses individual pileups at a specific position in the reference sequence.

        Args:
            position (int): The position in the reference sequence being analyzed.
            pileup: A pileup object representing a single base alignment from a read in the BAM file.

        The method calls 'process_snv' to handle potential SNVs at the given position. It is structured to be extendable for handling INDELs in the future.

        Note:
            Currently, the method for processing INDELs is commented out and not in use.

        Example:
            >>> # Assuming 'caller' is an instance of LiveVariantCaller and 'pileup' is a pileup object
            >>> caller.process_pileup_at_position(12345, pileup)
        """
        self.process_snv(position, pileup)
        # self.process_indel(position, pileup)

    def process_snv(self, position, pileup):
        """
        Processes a pileup for single nucleotide variant (SNV) detection.

        This method analyzes the pileup to identify potential SNVs at the specified position.

        Args:
            position (int): The position in the reference sequence where the variant is being checked.
            pileup: The pileup object from the BAM file representing a single read alignment at this position.

        If the pileup does not represent a deletion or a reference skip, it checks for SNVs. The method updates the 'memory' attribute with the detected SNV and its corresponding base quality score.

        Note:
            The method only considers non-deletion and non-reference skip pileups for SNV detection.

        Example:
            >>> # Assuming 'caller' is an instance of LiveVariantCaller and 'pileup' is a pileup object
            >>> caller.process_snv(12345, pileup)
        """
        if not pileup.is_del and not pileup.is_refskip:
            snv = pileup.alignment.query_sequence[pileup.query_position]

            if snv not in self.memory[position][c.VCF_SNVS].keys():
                self.memory[position][c.VCF_SNVS][snv] = []

            self.memory[position][c.VCF_SNVS][snv].append(pileup.alignment.query_qualities[pileup.query_position])

    def process_indel(self, position, pileup):
        """
        Processes a pileup for insertion and deletion (INDEL) detection.

        This method is designed to identify potential INDELs at the specified position in the reference sequence.

        Args:
            position (int): The position in the reference sequence where the INDEL is being checked.
            pileup: The pileup object from the BAM file representing a single read alignment at this position.

        The method checks if the pileup represents a deletion or a reference skip. If so, it updates the 'memory' attribute with the detected INDEL and, if available, its corresponding base quality score.

        Note:
            This method currently only handles deletions and reference skips. Insertions are marked for future implementation.

        Example:
            >>> # Assuming 'caller' is an instance of LiveVariantCaller and 'pileup' is a pileup object
            >>> caller.process_indel(12345, pileup)
        """
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
        """
        Processes the collected data in memory to prepare a list of identified variants.

        This method goes through the variant information stored in the 'memory' attribute and applies various filters and calculations to finalize the list of variants. It handles both SNVs (Single Nucleotide Variants) and INDELs (Insertions and Deletions).

        The method calculates statistics such as genotype likelihoods and quality scores for each potential variant. It then applies constraints based on minimum allele depth and evidence ratio to determine the final list of variants.

        Returns:
            List[dict]: A list of dictionaries, each representing a variant. The dictionary includes details such as variant position, alleles, quality scores, and other relevant information.

        The method uses a progress bar (tqdm) to display the progress of processing.

        Note:
            This method is the final step in the variant calling process and is expected to be called after all relevant BAM files have been processed and their data stored in memory.

        Example:
            >>> caller = LiveVariantCaller("reference.fasta", 20, 30, 10, 5, 0.1, 1000)
            >>> # After processing BAM files
            >>> variants = caller.prepare_variants()

        This method centralizes the statistical analysis and filtering logic, abstracting these complexities away from the main variant calling process.
        """
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
                            gl = genotype_likelihood
                            pl = u.to_phred_score(gl) #round(-10.0 * gl)
                        else:
                            gl = 0
                            pl = 0

                        # MAGIC HAPPENS HERE TOO
                        score = u.to_phred_score(1.0 - (genotype_likelihoods[allele] / sum_genotype_likelihoods))
                        
                        #qual = np.round(np.mean(snvs[allele])*100000, 2)

                        qual = u.to_phred_score(genotype_likelihoods[allele]) 
                        
                        if qual > 99:
                            qual = 64
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
        """
        Finds the previous variant in the list that is adjacent to the given variant.

        Args:
            variants (List[Variant]): A list of variants.
            variant (Variant): The variant for which the previous adjacent variant is being sought.

        Returns:
            Variant: The adjacent previous variant if found, otherwise None.

        This method is used for identifying a variant that is directly before the given variant in the genomic sequence. It's particularly useful in processes that require analysis of consecutive variants.

        Example:
            >>> adjacent_prev = caller.prev_variant(variants_list, current_variant)
        """
        return next(
            (
                v for v in variants
                if v[c.VCF_START] == variant[c.VCF_START] - 1
            ),
            None
        )

    def next_variant(self, variants: List[Variant], variant: Variant):
        """
        Finds the next variant in the list that is adjacent to the given variant.

        Args:
            variants (List[Variant]): A list of variants.
            variant (Variant): The variant for which the next adjacent variant is being sought.

        Returns:
            Variant: The adjacent next variant if found, otherwise None.

        This method is useful for identifying a variant that is directly after the given variant in the genomic sequence. It aids in processes where consecutive variant analysis is required.

        Example:
            >>> adjacent_next = caller.next_variant(variants_list, current_variant)
        """
        return next(
            (
                v for v in variants
                if v[c.VCF_START] == variant[c.VCF_START] + 1
            ),
            None
        )

    def concat_deletions(self, variants: List[Variant]):
        """
        Concatenates consecutive deletion variants in the list of variants.

        Args:
            variants (List[Variant]): A list of variants, potentially including consecutive deletions.

        Returns:
            List[Variant]: A new list of variants where consecutive deletions have been merged into single deletion events.

        This method is important for handling scenarios where multiple adjacent deletion variants may actually represent a single, larger deletion event. It iterates through the variants, merging adjacent deletions as needed.

        Example:
            >>> merged_variants = caller.concat_deletions(variants_list)
        """
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
        """
        Currently a placeholder for concatenating insertion variants.

        Args:
            variants (List[Variant]): A list of variants, which might include insertions.

        Returns:
            List[Variant]: Currently, this method returns the input list of variants unmodified.

        This method is intended to handle the concatenation of consecutive insertion variants. However, at this stage, it is a placeholder and does not perform any operations on the input list.

        Example:
            >>> merged_variants = caller.concat_insertions(variants_list)
        """
        return variants
    
    def write_vcf(self, output_vcf: str):
            
            """
        Writes the identified variants to a VCF (Variant Call Format) file.

        This method outputs the list of variants processed by the `LiveVariantCaller` into a VCF file, which is a standard format for representing variant calls in bioinformatics.

        Args:
            output_vcf (str): The file path where the VCF file will be written.

        The method initializes a VCF header using the pysam library, adding necessary metadata for each type of information included in the VCF file. The header includes metadata for total depth, allele depth, genotype likelihoods, phred-scaled genotype likelihoods, and a custom scoring function.

        For each reference sequence in the FASTA file, the method adds a contig entry to the VCF header.

        It then processes and writes each variant to the VCF file, sorting them by their start position and custom score.

        Note:
            The method assumes that variant calling and processing have been completed, and the variants are ready to be written to a file.

        Example:
            >>> caller = LiveVariantCaller("reference.fasta", 20, 30, 10, 5, 0.1, 1000)
            >>> caller.write_vcf("output.vcf")

        This method is crucial for generating a standard output that can be used for further analysis or reporting in genomic studies.
            """
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

    def create_checkpoint(self, filename):
        log.print_and_log(f'Creating checkpoint {filename}', log.INFO)
        print('SELF.MEMORY', type(self.memory))
        file = open(filename, 'wb')
        pickle.dump(self.memory, file)
        file.close()

    def load_checkpoint(self, filename):
        log.print_and_log(f'Loading checkpoint {filename}', log.INFO)

        file = open(filename, 'rb')
        self.memory = pickle.load(file)
        file.close()
    