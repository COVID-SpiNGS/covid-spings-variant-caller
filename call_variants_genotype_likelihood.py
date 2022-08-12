import pysam
import numpy as np
from config import minDepth
from utils import genotype_likelihood, is_relevant_position, to_error_probability


fastaFile = pysam.FastaFile('input/reference.fasta')
bamFile = pysam.AlignmentFile('input/input.bam', 'rb')
fileContent = ''

# should be 30175
pileupColumns = bamFile.pileup(min_base_quality=13)

fileContent += '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
    'POS', 
    'REF', 
    'NUM', 
    'A', 'A%',
    'T', 'T%',
    'C', 'C%',
    'G', 'G%'
)


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

        if(is_relevant_position(positions[-1])):
            fileContent += '%i\t%s\t%i\t%i\t%f\t%i\t%f\t%i\t%f\t%i\t%f\n' % (
                positions[-1]['pos'] + 1, 
                positions[-1]['ref'], 
                positions[-1]['numReads'], 
                len(positions[-1]['errorProbabilities']['A']), genotype_likelihood('A', positions[-1]),
                len(positions[-1]['errorProbabilities']['T']), genotype_likelihood('T', positions[-1]),
                len(positions[-1]['errorProbabilities']['C']), genotype_likelihood('C', positions[-1]),
                len(positions[-1]['errorProbabilities']['G']), genotype_likelihood('G', positions[-1]),
            )

bamFile.close()
fastaFile.close()

print(fileContent)

