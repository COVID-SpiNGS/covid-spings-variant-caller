
import pysam
import numpy as np

minDepth = 10
fastaFile = pysam.FastaFile('input/reference.fasta')
bamFile = pysam.AlignmentFile('input/input.bam', 'rb')
vcfHeader = pysam.VariantHeader()

for reference in fastaFile.references:
    vcfHeader.contigs.add(
        reference,
        fastaFile.get_reference_length(reference)
    )
    
vcfFile = pysam.VariantFile('output/call_variants_naive_vcf.vcf', mode='w', header=vcfHeader)
pileupColumns = bamFile.pileup()

def get_alternative(position, minNumReads=minDepth):
    if(position['numReads'] >= minNumReads):
        alt = max(position['calls'], key=position['calls'].get)
        errorProbability = 1.0 - (position['calls'][alt] / float(position['numReads']))
        phredQualityScore = -10  * np.log10(errorProbability) if errorProbability > 0.0 else 100

        return {
            'isRelevant': alt != position['ref'],
            'alt': alt,
            'qual': phredQualityScore
        }
    else:
        return {
            'isRelevant': False,
            'alt': None,
            'qual': 0.0
        }

positions = []
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
            }
        })

        for pileup in pileupColumn.pileups:
            if not pileup.is_del and not pileup.is_refskip:
                alt = pileup.alignment.query_sequence[pileup.query_position]
                positions[-1]['calls'][alt] += 1
                positions[-1]['numReads'] += 1


        alternative = get_alternative(positions[-1])
        if alternative['isRelevant']:
            vcfRecord = vcfFile.new_record(
                start=positions[-1]['pos'], 
                stop=positions[-1]['pos'] + 1,
                alleles=(
                    positions[-1]['base'], 
                    alternative['alt']
                ),
                qual=alternative['qual']
            )

            vcfFile.write(vcfRecord)
            
bamFile.close()
fastaFile.close()
vcfFile.close()