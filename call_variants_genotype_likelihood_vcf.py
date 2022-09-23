
from typing import List

import pysam
import functools
import operator
import numpy as np


from config import minRawDepth, minAlleleDepth, minMappingQuality, minBaseQuality
from structs import Position
from utils import genotype_likelihood, error_probability, to_error_probability, to_genotype_quality



fastaFile = pysam.FastaFile('input/reference.fasta')
bamFile = pysam.AlignmentFile('input/input.bam', 'rb')
vcfHeader = pysam.VariantHeader()

vcfHeader.add_meta('INFO', items=[
    ('ID', 'DP'),
    ('Number', 1),
    ('Type', 'Integer'),
    ('Description', 'Raw Depth')
])

vcfHeader.add_meta('INFO', items=[
    ('ID', 'AD'),
    ('Number', 1),
    ('Type', 'Integer'),
    ('Description', 'Allele Depth')
])

vcfHeader.add_meta('INFO', items=[
    ('ID', 'MBQ'),
    ('Number', 4),
    ('Type', 'Float'),
    ('Description', 'Mean Allele Base Qualities (A, T, C, G)')
])

vcfHeader.add_meta('INFO', items=[
    ('ID', 'MMQ'),
    ('Number', 4),
    ('Type', 'Float'),
    ('Description', 'Mean Allele Mapping Qualities (A, T, C, G)')
])

vcfHeader.add_meta('INFO', items=[
    ('ID', 'AF'),
    ('Number', 4),
    ('Type', 'Integer'),
    ('Description', 'Allele Frequency (A, T, C, G)')
])

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

def get_alternative(position, minAlleleDepth=minAlleleDepth):
    if(position['alleleDepth'] >= minAlleleDepth):
        alt = max(position['alleleFrequency'], key=position['alleleFrequency'].get)
        # qual = to_genotype_quality(genotype_likelihood(alt, position, 'alleleErrorProbabilities'))  
        qual = error_probability(alt, position, 'alleleErrorProbabilities')

        # if position['pos'] == 18873:
            # print(position)
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


        # print('genotypeLikelihoods', genotypeLikelihoods)
        if position['pos'] == 18873:
            print("")
           # print("")
            #print("")
            # print('genotypeLikelihoods', genotypeLikelihoods)
            # print('genotypeLikelihoodRatios', genotypeLikelihoodRatios)
            # print('pValues', pValues)


           # print("")
           # print("")
            print("")

            # print('pValues', pValues)
        # print('pValues', pValues)

        # alt = min(pValues, key=pValues.get)
        # qual = pValues[alt]
        # qual = to_genotype_quality(pValues[alt])  

        # How to go from likelihoodRatio to p-value
            

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

    if rawDepth >= minRawDepth:
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
                # errorProbability = 0.01

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


def sigmoid(x):
    # return 0.5
    return 1.0 / (1.0 + np.exp(-x))

def standard_score(x, mean, std):
    if(std > 0):
        return (x - mean) / std
    else: 
        return x

def feature_scaling(x, a=0.001, b=0.999):
    xMin = np.min(x)
    xMax = np.max(x)

    return (b - a) * (x - xMin) / (xMax - xMin) + a



for position in positions:

    allAlleleErrorProbabilities = np.concatenate((
        position['alleleErrorProbabilities']['A'], 
        position['alleleErrorProbabilities']['T'],
        position['alleleErrorProbabilities']['C'],
        position['alleleErrorProbabilities']['G']
    ), axis=0)

    allMean = np.mean(allAlleleErrorProbabilities) if len(allAlleleErrorProbabilities) else 0.0
    allStd = np.std(allAlleleErrorProbabilities) if len(allAlleleErrorProbabilities) else 1.0

    ## print(allAlleleErrorProbabilities)
    ## print(allMean, allStd)

    
    position['recalibratedAlleleErrorProbabilities'] = {
        'A': sigmoid(standard_score(position['alleleErrorProbabilities']['A'], allMean, allStd)),
        'C': sigmoid(standard_score(position['alleleErrorProbabilities']['C'], allMean, allStd)),
        'T': sigmoid(standard_score(position['alleleErrorProbabilities']['T'], allMean, allStd)),
        'G': sigmoid(standard_score(position['alleleErrorProbabilities']['G'], allMean, allStd)),
    }

    # print(position['alleleErrorProbabilities'])
    # print(position['recalibratedAlleleErrorProbabilities'



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