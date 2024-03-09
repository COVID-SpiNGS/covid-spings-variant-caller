import functools
import operator
import pickle
from time import localtime, strftime
from typing import List

import pysam
from pysam import AlignedSegment
from tqdm import tqdm

import src.variant_caller.utils as u
import src.variant_caller.vcf_file_constants as c
from src.config_util import logging as log
from src.variant_caller.models import Site, Variant


class LiveVariantCaller:
    """
    This class provides functionalities to identify variants based on various criteria such as base quality,
    mapping quality, allele depth, and evidence ratio. It uses a reference FASTA file for comparison
    and maintains an in-memory record of variant sites.

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

    def __init__(
        self,
        reference_fasta: str,
        min_base_quality: int,
        min_mapping_quality: int,
        min_total_depth: int,
        min_allele_depth: int,
        min_evidence_ratio: float,
        max_variants: int,
    ):
        """
        Constructor that sets up the LiveVariantCaller with the necessary parameters and resources for variant calling.

        @param reference_fasta: The file path to the reference genome FASTA file. This file is used to
        compare against sequencing data for variant calling.
        @type reference_fasta: str
        @param min_base_quality: The minimum quality score required for a base to be considered reliable.
        Quality scores lower than this threshold are ignored.
        @type min_base_quality: int
        @param min_mapping_quality: The minimum mapping quality score required for a read.
        Reads with mapping quality below this threshold are disregarded in variant calling.
        @type min_mapping_quality: int
        @param min_total_depth: The minimum number of reads ( = depth) that must cover a site for the site to be considered.
        @type min_total_depth: int
        @param min_allele_depth: The minimum number of reads supporting an allele to consider it a variant.
        @type min_allele_depth: int
        @param min_evidence_ratio: The minimum ratio of the number of reads supporting an allele to the total number
        of reads covering a site, required to consider it a variant.
        @type min_evidence_ratio: float
        @param max_variants: The maximum number of variant sites to keep in memory.
        @type max_variants: int
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

        The method closes the FASTA file associated with the instance, memory leaks and that the file is properly
        closed when the LiveVariantCaller object is no longer in use.

        """
        self.fasta_file.close()

    def reset_memory(self):
        """
        Resets the memory of the LiveVariantCaller instance.

        This method clears the current memory used for tracking variant sites. It initializes the 'memory' attribute to
        an empty dictionary.Each entry in this dictionary represents a variant site.

        """
        self.memory: dict[int, Site] = {}

    def process_bam(self, input_bam: str, reference_index=0):
        """
        Processes a BAM file to identify variant sites.

        This method reads a BAM file, iterates over its pileup columns, and processes each column to identify potential
        variant sites. The method utilizes the configured thresholds for base and mapping quality to filter reads.

        @param input_bam: File path to the input BAM file. This file contains the aligned sequencing reads.
        @type input_bam: str
        @param reference_index: The index of the reference sequence in the reference FASTA file to be used for pileup.
        Defaults to 0, this parameter is optional.
        @type reference_index: int, optional
        """

        bam_file = pysam.AlignmentFile(input_bam, "rb")
        pileup_columns = bam_file.pileup(
            min_mapping_quality=self.min_mapping_quality,
            min_base_quality=self.min_base_quality,
            reference=self.fasta_file.references[reference_index],
        )

        timestamp = strftime("[%Y-%m-%d %H:%M:%S]", localtime())
        progress_bar = tqdm(
            pileup_columns,
            desc=f"{timestamp} Processing {input_bam}",
            total=bam_file.get_reference_length(
                self.fasta_file.references[reference_index]
            ),
        )

        for pileup_column in progress_bar:
            self.process_pileup_column(pileup_column)

        bam_file.close()

    def process_pileup_column(self, pileup_column: AlignedSegment):
        """
        Processes a single pileup column from a BAM file to identify and record variant information.

        This method examines each pileup column in a BAM file to evaluate the presence of potential variants
        (SNVs and INDELs) at that specific position in the reference sequence.
        It calculates the total depth at the column's position and updates or initializes a memory structure
        with variant information for each position.

        @param pileup_column: A pileup column object representing a single position in the reference sequence, along
        with all the reads aligned to that position.
        @type pileup_column: pysam.AlignedSegment
        """

        total_depth = len(pileup_column.pileups)

        if pileup_column.reference_pos not in self.memory:
            reference = self.fasta_file.fetch(reference=pileup_column.reference_name)

            self.memory[pileup_column.reference_pos] = {
                c.VCF_REFERENCE: reference[pileup_column.reference_pos],
                c.VCF_TOTAL_DEPTH_KEY: total_depth,
                c.VCF_SNVS: {},
                c.VCF_INDELS: {},
            }
        else:
            self.memory[pileup_column.reference_pos][c.VCF_TOTAL_DEPTH_KEY] += (
                total_depth
            )

        for pileup in pileup_column.pileups:
            self.process_pileup_at_position(pileup_column.reference_pos, pileup)

    def process_pileup_at_position(self, position: int, pileup):
        """
        Processes a single pileup at a given position to identify single nucleotide variants (SNVs).

        @param position: The position in the reference sequence being analyzed. This is used to locate the specific base in the reference against which the reads are compared.
        @type position: int
        @param pileup: A pileup object representing a single base alignment from a read in the BAM file. This object encapsulates the alignment information for a read at the specified position, enabling the analysis of base mismatches that may indicate a variant.
        @type pileup: pysam.PileupRead or a similar custom object tailored to represent pileup information
        """
        self.process_snv(position, pileup)
        # self.process_indel(position, pileup)

    def process_snv(self, position, pileup):
        """
        Analyzes a pileup to identify potential single nucleotide variants (SNVs) at a specified position.

        @param position: The position in the reference sequence being analyzed for variants. This integer represents a specific base in the reference genome, and the analysis focuses on detecting variations at this exact location.
        @type position: int
        @param pileup: A pileup object representing the alignment of a single read against the reference sequence at the specified position. This object is essential for the analysis as it contains the alignment information needed to identify mismatches between the read and the reference sequence.
        @type pileup: pysam.PileupColumn or a similar object
        """
        if not pileup.is_del and not pileup.is_refskip:
            snv = pileup.alignment.query_sequence[pileup.query_position]

            if snv not in self.memory[position][c.VCF_SNVS].keys():
                self.memory[position][c.VCF_SNVS][snv] = []

            self.memory[position][c.VCF_SNVS][snv].append(
                pileup.alignment.query_qualities[pileup.query_position]
            )

    def process_indel(self, position, pileup):
        """
        Identifies potential insertions and deletions (INDELs) at a specified position in the reference sequence.

        This method examines the pileup at a given position for characteristics indicative of INDELs,
        such as deletions or reference skips.

        @param position: The position in the reference sequence being analyzed for INDELs.
        @type position: int
        @param pileup: A pileup object from the BAM file, representing the alignment of a single read against the reference sequence at the specified position. This object provides the necessary alignment details to discern the presence of INDELs.
        @type pileup: pysam.PileupColumn or a similar object that encapsulates the concept of a sequence aligment pileup.
        """
        if pileup.is_del or pileup.is_refskip:
            indel = (
                "-"
                if pileup.is_del
                else f"+{pileup.alignment.query_sequence[pileup.query_position]}"
            )

            if indel not in self.memory[position][c.VCF_INDELS].keys():
                # We could store more information here in the memory. But as the Base Qualty is the the only information that matters for further calculation we
                # save memory and only store them

                self.memory[position][c.VCF_INDELS][indel] = []

            if pileup.is_refskip:
                self.memory[position][c.VCF_INDELS][indel].append(
                    pileup.alignment.query_qualities[pileup.query_position]
                )
            else:
                self.memory[position][c.VCF_INDELS][indel].append(None)

    def prepare_variants(self):
        """
        Processes stored variant information to prepare a final list of identified variants.

        This method iterates through the 'memory' attribute, which contains preliminary variant data, and
        pplies a series of filters and statistical calculations. The aim is to refine this data into a finalized list of variants,
        encompassing both SNVs (Single Nucleotide Variants) and INDELs (Insertions and Deletions).

        @return: A list of dictionaries, with each dictionary representing a finalized variant. Each dictionary contains detailed information about the variant, including its position, alleles, quality scores, and additional relevant data.
        @rtype: List[dict]
        """
        timestamp = strftime("[%Y-%m-%d %H:%M:%S]", localtime())
        progress_bar = tqdm(
            self.memory,
            desc=f"{timestamp} Calculating statistics",
            total=len(self.memory.keys()),
        )

        variants: List[Variant] = []

        for position in progress_bar:
            if self.memory[position][c.VCF_TOTAL_DEPTH_KEY] >= self.min_total_depth:
                snvs = {
                    allele: [
                        quality for quality in self.memory[position][c.VCF_SNVS][allele]
                    ]
                    for allele in self.memory[position][c.VCF_SNVS].keys()
                }

                snvs_tuples = [(key, value) for key in snvs for value in snvs[key]]

                # MAGIC HAPPENS HERE
                if len(snvs_tuples) > 0:
                    x = u.get_likelihood(snvs_tuples)

                    genotype_likelihoods = u.extract_base_likelihood(
                        x, snvs_tuples, snvs
                    )

                    sum_genotype_likelihoods = functools.reduce(
                        operator.add, genotype_likelihoods.values(), 0.0
                    )
                    sum_genotype_likelihoods = (
                        sum_genotype_likelihoods
                        if sum_genotype_likelihoods != 0
                        else 1.0
                    )

                for allele in snvs.keys():
                    allele_depth = len(snvs[allele])

                    filter_constrains = [
                        self.memory[position][c.VCF_REFERENCE] != allele,
                        allele_depth >= self.min_allele_depth,
                        allele_depth / self.memory[position][c.VCF_TOTAL_DEPTH_KEY]
                        >= self.min_evidence_ratio,
                    ]

                    if all(filter_constrains):
                        genotype_likelihood = genotype_likelihoods[allele]

                        if genotype_likelihood != 0:
                            gl = genotype_likelihood
                            pl = u.to_phred_score(gl)  # round(-10.0 * gl)
                        else:
                            gl = 0
                            pl = 0

                        # MAGIC HAPPENS HERE TOO - But no idea how and why this happens?
                        score = u.to_phred_score(
                            1.0
                            - (genotype_likelihoods[allele] / sum_genotype_likelihoods)
                        )

                        # Very basic implementation
                        qual = u.to_phred_score(genotype_likelihoods[allele])

                        variants.append(
                            {
                                c.VCF_START: position,
                                c.VCF_STOP: position + 1,
                                c.VCF_ALLELES: (
                                    self.memory[position][c.VCF_REFERENCE],
                                    allele,
                                ),
                                c.VCF_QUAL: qual,
                                c.VCF_INFO: {
                                    c.VCF_DP: self.memory[position][
                                        c.VCF_TOTAL_DEPTH_KEY
                                    ],
                                    c.VCF_AD: allele_depth,
                                    c.VCF_GL: gl,
                                    c.VCF_PL: pl,
                                    c.VCF_SCORE: score,
                                },
                            }
                        )

                for indel in self.memory[position][c.VCF_INDELS].keys():
                    allele_depth = len(self.memory[position][c.VCF_INDELS][indel])

                    filter_constrains = [
                        allele_depth >= self.min_allele_depth,
                        allele_depth / self.memory[position][c.VCF_TOTAL_DEPTH_KEY]
                        >= self.min_evidence_ratio,
                    ]

                    if all(filter_constrains):
                        if indel == "-":
                            variants.append(
                                {
                                    c.VCF_START: position,
                                    c.VCF_STOP: position + 1,
                                    c.VCF_ALLELES: (
                                        self.memory[position][c.VCF_REFERENCE],
                                        "*",
                                    ),
                                    c.VCF_QUAL: 0,
                                    c.VCF_INFO: {
                                        c.VCF_DP: self.memory[position][
                                            c.VCF_TOTAL_DEPTH_KEY
                                        ],
                                        c.VCF_AD: allele_depth,
                                        c.VCF_GL: 0,
                                        c.VCF_PL: 0,
                                        c.VCF_SCORE: 0,
                                    },
                                }
                            )
                        else:
                            variants.append(
                                {
                                    c.VCF_START: position,
                                    c.VCF_STOP: position + 1,
                                    c.VCF_ALLELES: ("*", indel[1:]),
                                    c.VCF_QUAL: 0,
                                    c.VCF_INFO: {
                                        c.VCF_DP: self.memory[position][
                                            c.VCF_TOTAL_DEPTH_KEY
                                        ],
                                        c.VCF_ED: allele_depth,
                                        c.VCF_GL: 0,
                                        c.VCF_PL: 0,
                                        c.VCF_SCORE: 0,
                                    },
                                }
                            )

        return variants

    def prev_variant(self, variants: List[Variant], variant: Variant):
        """
        This method scans through a provided list of variants to find the variant that immediately precedes
        and is adjacent to a given variant.

        @param variants: A list of variants through which the method searches to find the adjacent previous variant.
        @type variants: List[Variant]
        @param variant: The variant for which the method seeks the adjacent previous variant.
        @type variant: Variant

        @return: The variant that is immediately preceding and adjacent to the given variant, if such a variant exists in the list; otherwise, None.
        @rtype: Variant or None


        """
        return next(
            (v for v in variants if v[c.VCF_START] == variant[c.VCF_START] - 1), None
        )

    def next_variant(self, variants: List[Variant], variant: Variant):
        """
        This method searches through a list of variants to identify the variant that immediately follows
        and is adjacent to a given variant.

        @param variants: The list of variants to be searched for the next adjacent variant.
        @type variants: List[Variant]
        @param variant: The reference variant for which the search for the next adjacent variant is conducted.
        @type variant: Variant

        @return: Returns the next variant that is adjacent to the specified variant if such a variant exists; otherwise, returns None.
        @rtype: Variant or None
        """
        return next(
            (v for v in variants if v[c.VCF_START] == variant[c.VCF_START] + 1), None
        )

    def concat_deletions(self, variants: List[Variant]):
        """
        Merges consecutive deletion variants into single deletion events within a list of variants.

        @param variants: The list of variants to be processed, which may contain sequences of consecutive deletions.
        @type variants: List[Variant]

        @return: A list of variants, where sequences of consecutive deletions have been merged into single deletion
        events, preserving the order and context of non-deletion variants.
        @rtype: List[Variant]
        """
        concatinated_variants: List[Variant] = []
        current_variant: Variant = None

        for variant in variants:
            if variant[c.VCF_ALLELES][1] == "*":
                next_variant = self.next_variant(variants, variant)

                if next_variant:
                    if not current_variant:
                        current_variant = variant
                    else:
                        current_variant = {
                            c.VCF_START: current_variant[c.VCF_START],
                            c.VCF_STOP: variant[c.VCF_STOP],
                            c.VCF_ALLELES: (
                                f"{current_variant[c.VCF_ALLELES][0]}{variant[c.VCF_ALLELES][0]}",
                                "*",
                            ),
                            c.VCF_QUAL: variant[c.VCF_QUAL],  # must be combined
                            c.VCF_INFO: variant[c.VCF_INFO],  # must be combined
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
        Placeholder method for future implementation of concatenating insertion variants.

        This method is intended to analyze a list of variants for sequences of consecutive insertions
        and combine them into single insertion events.

        @param variants: The list of variants to be processed, potentially including sequences of insertions.
        @type variants: List[Variant]

        @return: As of the current implementation, this method returns the input list of variants unaltered.
        @rtype: List[Variant]
        """
        return variants

    def write_vcf(self, output_vcf: str):
        """
        Exports the processed list of variants to a VCF file.

        Utilizes the pysam library to create a VCF file, starting with the initialization of a VCF header that includes
        essential metadata. The metadata contains information required for a comprehensive representation of variant calls,
        such as total depth, allele depth, genotype likelihoods, phred-scaled genotype likelihoods, and a custom scoring function.
        Additionally, for each reference sequence present in the associated FASTA file, a corresponding contig entry is added to the VCF header.

        Each variant, sorted by its start position and a custom score, is then processed and written to the VCF file.

        @param output_vcf: The filepath where the VCF file will be created and saved. This parameter specifies the destination for the output file, including the name of the file.
        @type output_vcf: str
        """
        print("VCF output", output_vcf)
        vcf_header = pysam.VariantHeader()

        vcf_header.add_meta(
            c.VCF_INFO,
            items=[
                (c.VCF_ID, c.VCF_DP),
                (c.VCF_NUMBER, 1),
                (c.VCF_TYPE, c.VCF_TYPE_INTEGER),
                (c.VCF_DESCRIPTION, c.VCF_TOTAL_DEPTH_STR),
            ],
        )

        vcf_header.add_meta(
            c.VCF_INFO,
            items=[
                (c.VCF_ID, c.VCF_AD),
                (c.VCF_NUMBER, 1),
                (c.VCF_TYPE, c.VCF_TYPE_INTEGER),
                (c.VCF_DESCRIPTION, c.VCF_ALLELE_DEPTH),
            ],
        )

        vcf_header.add_meta(
            c.VCF_INFO,
            items=[
                (c.VCF_ID, c.VCF_GL),
                (c.VCF_NUMBER, 1),
                (c.VCF_TYPE, c.VCF_TYPE_FLOAT),
                (
                    c.VCF_DESCRIPTION,
                    "Genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields",
                ),
            ],
        )

        vcf_header.add_meta(
            c.VCF_INFO,
            items=[
                (c.VCF_ID, c.VCF_PL),
                (c.VCF_NUMBER, 1),
                (c.VCF_TYPE, c.VCF_TYPE_INTEGER),
                (
                    c.VCF_DESCRIPTION,
                    "The phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field)",
                ),
            ],
        )

        vcf_header.add_meta(
            c.VCF_INFO,
            items=[
                (c.VCF_ID, c.VCF_SCORE),
                (c.VCF_NUMBER, 1),
                (c.VCF_TYPE, c.VCF_TYPE_FLOAT),
                (c.VCF_DESCRIPTION, "Custom scoring function"),
            ],
        )

        for reference in self.fasta_file.references:
            vcf_header.contigs.add(
                reference, self.fasta_file.get_reference_length(reference)
            )

        vcf_file = pysam.VariantFile(output_vcf, mode="w", header=vcf_header)

        variants = self.prepare_variants()
        # gvariants = self.concat_deletions(variants)

        for index, variant in enumerate(
            sorted(
                variants,
                key=lambda variant: (
                    variant[c.VCF_START],
                    variant[c.VCF_INFO][c.VCF_SCORE],
                ),
            )
        ):
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
        """
        Saves the current state to a file as a checkpoint.

        This method creates a checkpoint by serializing the current state of `self.memory` to a specified file.

        @param filename: The name of the file where the checkpoint will be saved.
        @type filename: str
        """

        log.print_and_log(f"Creating checkpoint {filename}", log.INFO)
        print("SELF.MEMORY", type(self.memory))
        file = open(filename, "wb")
        pickle.dump(self.memory, file)
        file.close()

    def load_checkpoint(self, filename):
        """
        Restores the state from a specified checkpoint file.

        @param filename: The name of the file from which the checkpoint will be loaded. The file must exist and contain
        a previously serialized state of `self.memory`.
        @type filename: str
        """
        log.print_and_log(f"Loading checkpoint {filename}", log.INFO)
        file = open(filename, "rb")
        self.memory = pickle.load(file)
        file.close()
