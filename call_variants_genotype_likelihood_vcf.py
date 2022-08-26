from typing import List
import pysam
import numpy as np
from config import minRawDepth, minAlleleDepth, minMappingQuality, minBaseQuality
from structs import Position
from utils import genotype_likelihood, to_error_probability, to_phred_quality_score, to_genotype_quality

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
    ('ID', 'ABQ'),
    ('Number', 4),
    ('Type', 'Float'),
    ('Description', 'Mean Allele Base Qualities (A, T, C, G)')
])

vcfHeader.add_meta('INFO', items=[
    ('ID', 'AMQ'),
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
    min_base_quality=minBaseQuality
)

def get_alternative(position, minAlleleDepth=minAlleleDepth):
    if(position['alleleDepth'] >= minAlleleDepth):
        alt = max(position['alleleFrequency'], key=position['alleleFrequency'].get)
        qual = to_genotype_quality(genotype_likelihood(alt, position))  

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
            'pos': pileupColumn.reference_pos,
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
                errorProbability = to_error_probability(pileup.alignment.query_qualities[pileup.query_position])

                positions[-1]['alleleFrequency'][alt] += 1

                positions[-1]['alleleErrorProbabilities'][alt] = np.append(
                    positions[-1]['alleleErrorProbabilities'][alt], 
                    errorProbability
                )

                positions[-1]['alleleBaseQualities'][alt] = np.append(
                    positions[-1]['alleleBaseQualities'][alt], 
                    pileup.alignment.query_qualities[pileup.query_position]
                )

                positions[-1]['alleleMappingQualities'][alt] = np.append(
                    positions[-1]['alleleMappingQualities'][alt], 
                    pileup.alignment.mapping_quality
                )

        positions[-1]['alleleDepth'] = positions[-1]['alleleFrequency']['A'] \
            + positions[-1]['alleleFrequency']['T'] \
            + positions[-1]['alleleFrequency']['C'] \
            + positions[-1]['alleleFrequency']['G']

        alternative = get_alternative(positions[-1])
        if alternative['isRelevant']:
            vcfRecord = vcfFile.new_record(
                start=positions[-1]['pos'], 
                stop=positions[-1]['pos'] + 1,
                alleles=(
                    positions[-1]['ref'], 
                    alternative['alt']
                ),
                qual=alternative['qual'],
                info={
                    'DP': positions[-1]['rawDepth'],
                    'AD': positions[-1]['alleleDepth'],
                    'ABQ': [
                        np.mean(positions[-1]['alleleBaseQualities']['A']) if len(positions[-1]['alleleBaseQualities']['A']) > 0 else 0,
                        np.mean(positions[-1]['alleleBaseQualities']['T']) if len(positions[-1]['alleleBaseQualities']['T']) > 0 else 0,
                        np.mean(positions[-1]['alleleBaseQualities']['C']) if len(positions[-1]['alleleBaseQualities']['C']) > 0 else 0,
                        np.mean(positions[-1]['alleleBaseQualities']['G']) if len(positions[-1]['alleleBaseQualities']['G']) > 0 else 0
                    ],
                    'AMQ': [
                        np.mean(positions[-1]['alleleMappingQualities']['A']) if len(positions[-1]['alleleMappingQualities']['A']) > 0 else 0,
                        np.mean(positions[-1]['alleleMappingQualities']['T']) if len(positions[-1]['alleleMappingQualities']['T']) > 0 else 0,
                        np.mean(positions[-1]['alleleMappingQualities']['C']) if len(positions[-1]['alleleMappingQualities']['C']) > 0 else 0,
                        np.mean(positions[-1]['alleleMappingQualities']['G']) if len(positions[-1]['alleleMappingQualities']['G']) > 0 else 0
                    ],
                    'AF': list(positions[-1]['alleleFrequency'].values())
                }
            )

            vcfFile.write(vcfRecord)

bamFile.close()
fastaFile.close()
vcfFile.close()