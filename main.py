from variant_caller.live_variant_caller import LiveVariantCaller
from os.path import getsize

import time


def main():
    # start = time.time()

    minEvidenceDepth = 5
    minEvidenceRatio = 0.0
    maxVariants = 1
    minTotalDepth = 10
    minMappingQuality = 20
    minBaseQuality = 30

    liveVariantCaller = LiveVariantCaller(
        'input/reference-covid.fasta',
        minBaseQuality,
        minMappingQuality,
        minTotalDepth,
        minEvidenceDepth,
        minEvidenceRatio,
        maxVariants
    )

    start = time.time()
    file = 'input/input-10/input.bam'
    liveVariantCaller.process_bam(file)
    print('Filesize', getsize(file))
    print('Time', time.time() - start)
    # liveVariantCaller.process_bam('input/input.bam')
    #
    # 
    # liveVariantCaller.process_bam('input/real/input_001.bam')
    # print('1', time.time() - start)
    # 
    # 
    # liveVariantCaller.process_bam('input/real/input_002.bam')
    # print('2', time.time() - start)
    # 
    # 
    # liveVariantCaller.process_bam('input/real/input_003.bam')
    # print('3', time.time() - start)
    # 
    # liveVariantCaller.process_bam('input/real/input_004.bam')
    # print('4', time.time() - start)
    # 
    # liveVariantCaller.process_bam('input/real/input_005.bam')
    # print('5', time.time() - start)
    # 
    # liveVariantCaller.process_bam('input/real/input_006.bam')
    # print('6', time.time() - start)
    # 
    # liveVariantCaller.process_bam('input/real/input_007.bam')
    # print('7', time.time() - start)
    # 
    # liveVariantCaller.process_bam('input/real/input_008.bam')
    # print('8', time.time() - start)
    # 
    # liveVariantCaller.process_bam('input/real/input_009.bam')
    # print('9', time.time() - start)
    # 
    # liveVariantCaller.process_bam('input/real/input_010.bam')
    # print('10', time.time() - start)
    # 
    # 
    # liveVariantCaller.write_vcf('output/live_multi_variant_caller.vcf')
    # # liveVariantCaller.create_checkpoint('output/checkpoint.pkl')
    #
    # print('total', time.time() - start)


if __name__ == '__main__':
    main()
