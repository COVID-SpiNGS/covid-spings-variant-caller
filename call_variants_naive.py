import pysam

fastaFile = pysam.FastaFile('input/reference.fasta')
bamFile = pysam.AlignmentFile('input/input.bam', 'rb')
min_depth = 10
probabilities = []
fileContent = ''

# should be 30175
pileupColumns = bamFile.pileup()
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
            }
        })

        for pileup in pileupColumn.pileups:
            if not pileup.is_del and not pileup.is_refskip:
                readBase = pileup.alignment.query_sequence[pileup.query_position]
                positions[-1]['calls'][readBase] += 1
                positions[-1]['num_reads'] += 1


        if is_relevant_position(positions[-1]):
            fileContent += '%i\t%s\t%i\t%i\t%f\t%i\t%f\t%i\t%f\t%i\t%f\n' % (
                positions[-1]['position'], 
                positions[-1]['base'], 
                positions[-1]['num_reads'],
                positions[-1]['calls']['A'], (positions[-1]['calls']['A'] / float(positions[-1]['num_reads'])),
                positions[-1]['calls']['T'], (positions[-1]['calls']['T'] / float(positions[-1]['num_reads'])),
                positions[-1]['calls']['C'], (positions[-1]['calls']['C'] / float(positions[-1]['num_reads'])),
                positions[-1]['calls']['G'], (positions[-1]['calls']['G'] / float(positions[-1]['num_reads']))
            )
            
bamFile.close()
fastaFile.close()

print(fileContent)

