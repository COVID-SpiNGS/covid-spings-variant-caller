import pysam
import numpy as np

bamFile = pysam.AlignmentFile('input/input.bam', 'rb')
min_depth = 10
probabilities = []
fileContent = ''

i = 0

# should be 30175
pileupColumns = bamFile.pileup(min_base_quality=13)

fileContent += '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
    'POS', 
    'REF', 
    'DEPTH', 
    'A', 'A%',
    'T', 'T%',
    'C', 'C%',
    'G', 'G%'
)

def quality_to_error_probability(quality: int):
    return np.power(10, quality / -10)


def likelyhood(base: str, position): 
    if base == 'A':
        return (1.0 - position['errorProbabilities']['A']).prod() \
            * position['errorProbabilities']['T'].prod() \
            * position['errorProbabilities']['C'].prod() \
            * position['errorProbabilities']['G'].prod()
    elif base == 'T':
        return (1.0 - position['errorProbabilities']['T']).prod() \
            * position['errorProbabilities']['A'].prod() \
            * position['errorProbabilities']['C'].prod() \
            * position['errorProbabilities']['G'].prod()
    elif base == 'C':
        return (1.0 - position['errorProbabilities']['C']).prod() \
            * position['errorProbabilities']['A'].prod() \
            * position['errorProbabilities']['T'].prod() \
            * position['errorProbabilities']['G'].prod()
    elif base == 'G':
        return (1.0 - position['errorProbabilities']['G']).prod() \
            * position['errorProbabilities']['A'].prod() \
            * position['errorProbabilities']['T'].prod() \
            * position['errorProbabilities']['C'].prod()


positions = []
for pileupColumn in pileupColumns:
    depth = len(pileupColumn.pileups)

    if depth >= min_depth:
        positions.append({
            'position': i,
            'base': '?',
            'depth': depth,
            'errorProbabilities': {
                'A': np.array([]),
                'T': np.array([]),
                'C': np.array([]),
                'G': np.array([])
            }
        })
        
        for pileup in pileupColumn.pileups:
            if not pileup.is_del and not pileup.is_refskip:
                readBase = pileup.alignment.query_sequence[pileup.query_position]
                errorProbability = quality_to_error_probability(pileup.alignment.query_qualities[pileup.query_position])

                positions[-1]['errorProbabilities'][readBase] = np.append(positions[-1]['errorProbabilities'][readBase], errorProbability)


        fileContent += '%i\t%s\t%i\t%i\t%f\t%i\t%f\t%i\t%f\t%i\t%f\n' % (
            positions[-1]['position'], 
            positions[-1]['base'], 
            positions[-1]['depth'], 
            len(positions[-1]['errorProbabilities']['A']), likelyhood('A', positions[-1]),
            len(positions[-1]['errorProbabilities']['T']), likelyhood('T', positions[-1]),
            len(positions[-1]['errorProbabilities']['C']), likelyhood('C', positions[-1]),
            len(positions[-1]['errorProbabilities']['G']), likelyhood('G', positions[-1]),
        )


    i = i + 1


print(fileContent)