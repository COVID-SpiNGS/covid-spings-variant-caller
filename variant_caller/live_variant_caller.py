import functools
import operator
import pysam
import numpy as np

from pysam import AlignedSegment
from structs import Alternative, Position

from utils import error_probability, genotype_likelihood, to_error_probability, to_phred_scale, to_phred_quality_score

class LiveVariantCaller:

    def __init__(self, referenceFasta: str, minBaseQuality: int, minMappingQuality: int, minTotalDepth: int, minCandidatesDepth: int):
        self.minBaseQuality = minBaseQuality
        self.minMappingQuality = minMappingQuality
        self.minTotalDepth = minTotalDepth
        self.minCandidatesDepth = minCandidatesDepth
        
        self.fastaFile = pysam.FastaFile(referenceFasta)

        self.reset_memory()

    def __del__(self):
        self.fastaFile.close()

    def process_bam(self, inputBam: str):
        bamFile = pysam.AlignmentFile(inputBam, 'rb')
        
        pileupColumns = bamFile.pileup(
            min_mapping_quality=self.minMappingQuality,
            min_base_quality=self.minBaseQuality,
        )

        for pileupColumn in pileupColumns:
            self.process_pileup_column(pileupColumn)

        bamFile.close()

    def process_pileup_column(self, pileupColumn: AlignedSegment):
        totalDepth = len(pileupColumn.pileups)

        if pileupColumn.reference_pos not in self.memory:
            reference = self.fastaFile.fetch(reference=pileupColumn.reference_name)

            self.memory[pileupColumn.reference_pos] = {
                'ref': reference[pileupColumn.reference_pos],
                'totalDepth': totalDepth,
                'baseFrequencies': {
                    'A': 0,
                    'T': 0,
                    'C': 0,
                    'G': 0,
                },
                'baseQualities': {
                    'A': np.array([]),
                    'T': np.array([]),
                    'C': np.array([]),
                    'G': np.array([])
                },
                'mappingQualities': {
                    'A': np.array([]),
                    'T': np.array([]),
                    'C': np.array([]),
                    'G': np.array([])
                },
                'baseErrorProbabilities': {
                    'A': np.array([]),
                    'T': np.array([]),
                    'C': np.array([]),
                    'G': np.array([])
                },
            }
        else:
            self.memory[pileupColumn.reference_pos]['totalDepth'] += totalDepth

        for pileup in pileupColumn.pileups:
            self.process_pileup_at_position(pileupColumn.reference_pos, pileup)

        self.memory[pileupColumn.reference_pos]['candidatesDepth'] = self.memory[pileupColumn.reference_pos]['baseFrequencies']['A'] \
            + self.memory[pileupColumn.reference_pos]['baseFrequencies']['T'] \
            + self.memory[pileupColumn.reference_pos]['baseFrequencies']['C'] \
            + self.memory[pileupColumn.reference_pos]['baseFrequencies']['G']
            

    def process_pileup_at_position(self, position: int, pileup):
        if not pileup.is_del and not pileup.is_refskip:
            alt = pileup.alignment.query_sequence[pileup.query_position]
            self.memory[position]['baseFrequencies'][alt] += 1


            self.memory[position]['baseQualities'][alt] = np.append(
                self.memory[position]['baseQualities'][alt], 
                pileup.alignment.query_qualities[pileup.query_position]
            )

            self.memory[position]['mappingQualities'][alt] = np.append(
                self.memory[position]['mappingQualities'][alt], 
                pileup.alignment.mapping_quality
            )


            self.memory[position]['baseErrorProbabilities'][alt] = np.append(
                self.memory[position]['baseErrorProbabilities'][alt], 
                to_error_probability(pileup.alignment.query_qualities[pileup.query_position])
            )


    def reset_memory(self):
        self.memory: dict[int, Position] = {}

    def get_alternative(self, position: Position) -> Alternative:
        if(position['candidatesDepth'] >= self.minCandidatesDepth):
            genotypeLikelihoods = {
                base: (
                    genotype_likelihood(base, position, 'baseErrorProbabilities') 
                    if len(position['baseErrorProbabilities'][base] > 0)
                    else 0.0
                )
                for base in position['baseErrorProbabilities'].keys()
            }

            sumGenotypeLikelihoods = functools.reduce(operator.add, genotypeLikelihoods.values())

            pValues = {
                base: (genotypeLikelihoods[base] / sumGenotypeLikelihoods)
                for base in genotypeLikelihoods.keys()
            }

            alt = max(position['baseFrequencies'], key=position['baseFrequencies'].get)
            qual = to_phred_scale(pValues[alt])  

            return {
                'isRelevant': alt != position['ref'],
                'alt': alt,
                'qual': qual
            }
        else:
            return {
                'isRelevant': False,
                'alt': None,
                'qual': 0
            }

    def write_vcf(self, outputVfc: str):
        vcfHeader = pysam.VariantHeader()

        vcfHeader.add_meta('INFO', items=[
            ('ID', 'TD'),
            ('Number', 1),
            ('Type', 'Integer'),
            ('Description', 'Total Depth')
        ])

        vcfHeader.add_meta('INFO', items=[
            ('ID', 'CD'),
            ('Number', 1),
            ('Type', 'Integer'),
            ('Description', 'Candidates Depth')
        ])

        vcfHeader.add_meta('INFO', items=[
            ('ID', 'MBQ'),
            ('Number', 4),
            ('Type', 'Float'),
            ('Description', 'Mean Base Qualities (A, T, C, G)')
        ])

        vcfHeader.add_meta('INFO', items=[
            ('ID', 'MMQ'),
            ('Number', 4),
            ('Type', 'Float'),
            ('Description', 'Mean Mapping Qualities (A, T, C, G)')
        ])

        vcfHeader.add_meta('INFO', items=[
            ('ID', 'BF'),
            ('Number', 4),
            ('Type', 'Integer'),
            ('Description', 'Base Frequencies (A, T, C, G)')
        ])

        for reference in self.fastaFile.references:
            vcfHeader.contigs.add(
                reference,
                self.fastaFile.get_reference_length(reference)
            )

        vcfFile = pysam.VariantFile(outputVfc, mode='w', header=vcfHeader)

        for index in self.memory:
            if self.memory[index]['totalDepth'] >= self.minTotalDepth and self.memory[index]['candidatesDepth'] >= self.minCandidatesDepth:
                alternative = self.get_alternative(self.memory[index])

                if alternative['isRelevant']:
                    vcfRecord = vcfFile.new_record(
                        start=index, 
                        stop=index + 1,
                        alleles=(
                            self.memory[index]['ref'], 
                            alternative['alt']
                        ),
                        qual=alternative['qual'],
                        info={
                            'TD': self.memory[index]['totalDepth'],
                            'CD': self.memory[index]['candidatesDepth'],
                            'MBQ': [
                                np.mean(self.memory[index]['baseQualities']['A']) if len(self.memory[index]['baseQualities']['A']) > 0 else 0,
                                np.mean(self.memory[index]['baseQualities']['T']) if len(self.memory[index]['baseQualities']['T']) > 0 else 0,
                                np.mean(self.memory[index]['baseQualities']['C']) if len(self.memory[index]['baseQualities']['C']) > 0 else 0,
                                np.mean(self.memory[index]['baseQualities']['G']) if len(self.memory[index]['baseQualities']['G']) > 0 else 0
                            ],
                            'MMQ': [
                                np.mean(self.memory[index]['mappingQualities']['A']) if len(self.memory[index]['mappingQualities']['A']) > 0 else 0,
                                np.mean(self.memory[index]['mappingQualities']['T']) if len(self.memory[index]['mappingQualities']['T']) > 0 else 0,
                                np.mean(self.memory[index]['mappingQualities']['C']) if len(self.memory[index]['mappingQualities']['C']) > 0 else 0,
                                np.mean(self.memory[index]['mappingQualities']['G']) if len(self.memory[index]['mappingQualities']['G']) > 0 else 0
                            ],
                            'BF': list(self.memory[index]['baseFrequencies'].values())
                        }
                    )

                    vcfFile.write(vcfRecord)
        vcfFile.close()
