
from typing import List
import pysam


bamFile = pysam.AlignmentFile('input.input.bam', 'rb')
min_depth = 10
probabilities = []
fileContent = ''

i = 0

# should be 30175
pileupColumns = bamFile.pileup()

fileContent += '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
    'POS', 
    'REF', 
    'DEPTH', 
    'A', 'A%',
    'T', 'T%',
    'C', 'C%',
    'G', 'G%'
)

positions = []
for pileupColumn in pileupColumns:
    depth = len(pileupColumn.pileups)


    if depth >= min_depth:
        positions.append({
            'position': i,
            'base': '?',
            'depth': depth,
            'calls': {
                'A': 0,
                'T': 0,
                'C': 0,
                'G': 0,
            }
        })
        
        for pileup in pileupColumn.pileups:
            referenceBase = '?'


            if not pileup.is_del and not pileup.is_refskip:
                readBase = pileup.alignment.query_sequence[pileup.query_position]
                positions[-1]['calls'][readBase] += 1


        totalCalls = positions[-1]['calls']['A'] + positions[-1]['calls']['T'] + positions[-1]['calls']['G'] + positions[-1]['calls']['C']


        fileContent += '%i\t%s\t%i\t%i\t%f\t%i\t%f\t%i\t%f\t%i\t%f\n' % (
            i, 
            referenceBase, 
            depth, 
            positions[-1]['calls']['A'], (positions[-1]['calls']['A'] / float(totalCalls)) if totalCalls > 0 else 0,
            positions[-1]['calls']['T'], (positions[-1]['calls']['T'] / float(totalCalls)) if totalCalls > 0 else 0,
            positions[-1]['calls']['C'], (positions[-1]['calls']['C'] / float(totalCalls)) if totalCalls > 0 else 0,
            positions[-1]['calls']['G'], (positions[-1]['calls']['G'] / float(totalCalls)) if totalCalls > 0 else 0
        )

    i = i + 1

print(fileContent)