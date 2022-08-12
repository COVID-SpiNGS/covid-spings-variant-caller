import pysam

minDepth = 10
fastaFile = pysam.FastaFile('input/reference.fasta')
bamFile = pysam.AlignmentFile('input/input.bam', 'rb')
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


def is_relevant_position(position, minNumReads=minDepth):
    if(position['numReads'] >= minNumReads):
        alt = max(position['calls'], key=position['calls'].get)
        return alt != position['base']
    else:
        return False

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


        if is_relevant_position(positions[-1]):
            fileContent += '%i\t%s\t%i\t%i\t%f\t%i\t%f\t%i\t%f\t%i\t%f\n' % (
                positions[-1]['pos'] + 1, 
                positions[-1]['ref'], 
                positions[-1]['numReads'],
                positions[-1]['calls']['A'], (positions[-1]['calls']['A'] / float(positions[-1]['numReads'])),
                positions[-1]['calls']['T'], (positions[-1]['calls']['T'] / float(positions[-1]['numReads'])),
                positions[-1]['calls']['C'], (positions[-1]['calls']['C'] / float(positions[-1]['numReads'])),
                positions[-1]['calls']['G'], (positions[-1]['calls']['G'] / float(positions[-1]['numReads']))
            )
            
bamFile.close()
fastaFile.close()

print(fileContent)

