from typing import List
import pysam
import numpy as np
from config import minDepth
from structs import Position
from utils import genotype_likelihood, to_error_probability, to_phred_quality_score, to_genotype_quality

fastaFile = pysam.FastaFile('input/reference.fasta')
bamFile = pysam.AlignmentFile('input/input.bam', 'rb')
vcfHeader = pysam.VariantHeader()

for reference in fastaFile.references:
    vcfHeader.contigs.add(
        reference,
        fastaFile.get_reference_length(reference)
    )
    
vcfFile = pysam.VariantFile('output/call_variants_genotype_likelihood_vcf.vcf', mode='w', header=vcfHeader)
pileupColumns = bamFile.pileup(min_base_quality=13)

def get_alternative(position, minNumReads=minDepth):
    if(position['numReads'] >= minNumReads):
        alt = max(position['calls'], key=position['calls'].get)
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
    depth = len(pileupColumn.pileups)

    if depth >= minDepth:
        reference = fastaFile.fetch(reference=pileupColumn.reference_name)
        positions.append({
            'pos': pileupColumn.reference_pos,
            'ref': reference[pileupColumn.reference_pos],
            'numReads': 0,
            'calls': {
                'A': 0,
                'T': 0,
                'C': 0,
                'G': 0,
            },
            'errorProbabilities': {
                'A': np.array([]),
                'T': np.array([]),
                'C': np.array([]),
                'G': np.array([])
            },
            'phredQualityScores': {
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

                positions[-1]['calls'][alt] += 1
                positions[-1]['numReads'] += 1
                positions[-1]['errorProbabilities'][alt] = np.append(positions[-1]['errorProbabilities'][alt], errorProbability)
                positions[-1]['phredQualityScores'][alt] = np.append(positions[-1]['phredQualityScores'][alt], pileup.alignment.query_qualities[pileup.query_position])

        alternative = get_alternative(positions[-1])
        if alternative['isRelevant']:
            vcfRecord = vcfFile.new_record(
                start=positions[-1]['pos'], 
                stop=positions[-1]['pos'] + 1,
                alleles=(
                    positions[-1]['ref'], 
                    alternative['alt']
                ),
                qual=alternative['qual']
            )

            vcfFile.write(vcfRecord)

bamFile.close()
fastaFile.close()
vcfFile.close()