import pysam
import numpy as np

fastaFile = pysam.FastaFile('input/reference.fasta')
bamFile = pysam.AlignmentFile('input/input.bam', 'rb')
min_depth = 10
probabilities = []
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

def is_relevant_position(position, min_num_reads=min_depth):
    if(position['num_reads'] >= min_num_reads):
        readBase = max(position['calls'], key=position['calls'].get)
        return readBase != position['base']
    else:
        return False

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
        reference = fastaFile.fetch(reference=pileupColumn.reference_name)
        positions.append({
            'position': pileupColumn.reference_pos,
            'base': reference[pileupColumn.reference_pos],
            'num_reads': 0,
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
            }
        })
        
        for pileup in pileupColumn.pileups:
            if not pileup.is_del and not pileup.is_refskip:
                readBase = pileup.alignment.query_sequence[pileup.query_position]
                errorProbability = quality_to_error_probability(pileup.alignment.query_qualities[pileup.query_position])

                positions[-1]['calls'][readBase] += 1
                positions[-1]['num_reads'] += 1
                positions[-1]['errorProbabilities'][readBase] = np.append(positions[-1]['errorProbabilities'][readBase], errorProbability)

        if(is_relevant_position(positions[-1])):
            fileContent += '%i\t%s\t%i\t%i\t%f\t%i\t%f\t%i\t%f\t%i\t%f\n' % (
                positions[-1]['position'], 
                positions[-1]['base'], 
                positions[-1]['num_reads'], 
                len(positions[-1]['errorProbabilities']['A']), likelyhood('A', positions[-1]),
                len(positions[-1]['errorProbabilities']['T']), likelyhood('T', positions[-1]),
                len(positions[-1]['errorProbabilities']['C']), likelyhood('C', positions[-1]),
                len(positions[-1]['errorProbabilities']['G']), likelyhood('G', positions[-1]),
            )

bamFile.close()
fastaFile.close()

print(fileContent)

