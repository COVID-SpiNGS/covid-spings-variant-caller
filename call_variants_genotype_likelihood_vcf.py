
from typing import List

import pysam
import functools
import operator
import numpy as np


from config import minTotalDepth, minCandidatesDepth, minMappingQuality, minBaseQuality
from structs import Position
from utils import genotype_likelihood, error_probability, to_error_probability, to_phred_scale



fastaFile = pysam.FastaFile('input/reference.fasta')
bamFile = pysam.AlignmentFile('input/input.bam', 'rb')

for reference in fastaFile.references:
    vcfHeader.contigs.add(
        reference,
        fastaFile.get_reference_length(reference)
    )
    
vcfFile = pysam.VariantFile('output/call_variants_genotype_likelihood_vcf.vcf', mode='w', header=vcfHeader)

pileupColumns = bamFile.pileup(
    min_mapping_quality=minMappingQuality,
    min_base_quality=minBaseQuality,
)

def get_alternative(position, minAlleleDepth=minCandidatesDepth):
    if(position['alleleDepth'] >= minAlleleDepth):
        genotypeLikelihoods = {
            base: (
                genotype_likelihood(base, position, 'alleleErrorProbabilities') 
                if len(position['alleleErrorProbabilities'][base] > 0)
                else 0.0
            )
            for base in position['alleleErrorProbabilities'].keys()
        }

        sumGenotypeLikelihoods = functools.reduce(operator.add, genotypeLikelihoods.values())

        pValues = {
            base: 1.0 - (genotypeLikelihoods[base] / sumGenotypeLikelihoods)
            for base in genotypeLikelihoods.keys()
        }

        alt = max(position['alleleFrequency'], key=position['alleleFrequency'].get)
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
            'qual': 0.0
        }

positions: List[Position] = []
for pileupColumn in pileupColumns:
    rawDepth = len(pileupColumn.pileups)

    if rawDepth >= minTotalDepth:
        reference = fastaFile.fetch(reference=pileupColumn.reference_name)
        positions.append({
            'pos': pileupColumn.reference_pos + 1,
            'ref': reference[pileupColumn.reference_pos],
            'rawDepth': rawDepth,
            'alleleDepth': 0,
            'alleleFrequency': {
                'A': 0,
                'T': 0,
                'C': 0,
                'G': 0,
            },
            'alleleErrorProbabilities': {
                'A': np.array([]),
                'T': np.array([]),
                'C': np.array([]),
                'G': np.array([])
            },
            'recalibratedAlleleErrorProbabilities': {
                'A': np.array([]),
                'T': np.array([]),
                'C': np.array([]),
                'G': np.array([])
            },
            'alleleBaseQualities': {
                'A': np.array([]),
                'T': np.array([]),
                'C': np.array([]),
                'G': np.array([])
            },
            'alleleMappingQualities': {
                'A': np.array([]),
                'T': np.array([]),
                'C': np.array([]),
                'G': np.array([])
            }
        })
        
        for pileup in pileupColumn.pileups:
            if not pileup.is_del and not pileup.is_refskip:
                alt = pileup.alignment.query_sequence[pileup.query_position]
                phredQuality = pileup.alignment.query_qualities[pileup.query_position]
                errorProbability = to_error_probability(phredQuality)

                positions[-1]['alleleFrequency'][alt] += 1

                positions[-1]['alleleErrorProbabilities'][alt] = np.append(
                    positions[-1]['alleleErrorProbabilities'][alt], 
                    errorProbability
                )

                positions[-1]['alleleBaseQualities'][alt] = np.append(
                    positions[-1]['alleleBaseQualities'][alt], 
                    phredQuality
                )

                positions[-1]['alleleMappingQualities'][alt] = np.append(
                    positions[-1]['alleleMappingQualities'][alt], 
                    pileup.alignment.mapping_quality
                )

        positions[-1]['alleleDepth'] = positions[-1]['alleleFrequency']['A'] \
            + positions[-1]['alleleFrequency']['T'] \
            + positions[-1]['alleleFrequency']['C'] \
            + positions[-1]['alleleFrequency']['G']



for position in positions:
    alternative = get_alternative(position)

    if alternative['isRelevant']:
        vcfRecord = vcfFile.new_record(
            start=position['pos'] - 1, 
            stop=position['pos'],
            alleles=(
                position['ref'], 
                alternative['alt']
            ),
            qual=alternative['qual'],
            info={
                'DP': position['rawDepth'],
                'AD': position['alleleDepth'],
                'MBQ': [
                    np.mean(position['alleleBaseQualities']['A']) if len(position['alleleBaseQualities']['A']) > 0 else 0,
                    np.mean(position['alleleBaseQualities']['T']) if len(position['alleleBaseQualities']['T']) > 0 else 0,
                    np.mean(position['alleleBaseQualities']['C']) if len(position['alleleBaseQualities']['C']) > 0 else 0,
                    np.mean(position['alleleBaseQualities']['G']) if len(position['alleleBaseQualities']['G']) > 0 else 0
                ],
                'MMQ': [
                    np.mean(position['alleleMappingQualities']['A']) if len(position['alleleMappingQualities']['A']) > 0 else 0,
                    np.mean(position['alleleMappingQualities']['T']) if len(position['alleleMappingQualities']['T']) > 0 else 0,
                    np.mean(position['alleleMappingQualities']['C']) if len(position['alleleMappingQualities']['C']) > 0 else 0,
                    np.mean(position['alleleMappingQualities']['G']) if len(position['alleleMappingQualities']['G']) > 0 else 0
                ],
                'AF': list(position['alleleFrequency'].values())
            }
        )

        vcfFile.write(vcfRecord)



bamFile.close()
fastaFile.close()
vcfFile.close()