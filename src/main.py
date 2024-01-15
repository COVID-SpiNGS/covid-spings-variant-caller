from live_variant_caller import LiveVariantCaller
from os.path import getsize

import time


def main():
    # start = time.time()

    min_evidence_depth = 5
    min_evidence_ratio = 0.0
    max_variants = 1
    min_total_depth = 10
    min_mapping_quality = 20
    min_base_quality = 30

    live_variant_caller = LiveVariantCaller(
        'input/reference-covid.fasta',
        min_base_quality,
        min_mapping_quality,
        min_total_depth,
        min_evidence_depth,
        min_evidence_ratio,
        max_variants
    )

    start = time.time()
    file = 'input/input.bam'
    live_variant_caller.process_bam(file)
    print('Filesize', getsize(file))
    print('Time', time.time() - start)
    # live_variant_caller.process_bam('input/input.bam')
    #
    # 
    # live_variant_caller.process_bam('input/real/input_001.bam')
    # print('1', time.time() - start)
    # 
    # 
    # live_variant_caller.process_bam('input/real/input_002.bam')
    # print('2', time.time() - start)
    # 
    # 
    # live_variant_caller.process_bam('input/real/input_003.bam')
    # print('3', time.time() - start)
    # 
    # live_variant_caller.process_bam('input/real/input_004.bam')
    # print('4', time.time() - start)
    # 
    # live_variant_caller.process_bam('input/real/input_005.bam')
    # print('5', time.time() - start)
    # 
    # live_variant_caller.process_bam('input/real/input_006.bam')
    # print('6', time.time() - start)
    # 
    # live_variant_caller.process_bam('input/real/input_007.bam')
    # print('7', time.time() - start)
    # 
    # live_variant_caller.process_bam('input/real/input_008.bam')
    # print('8', time.time() - start)
    # 
    # live_variant_caller.process_bam('input/real/input_009.bam')
    # print('9', time.time() - start)
    # 
    # live_variant_caller.process_bam('input/real/input_010.bam')
    # print('10', time.time() - start)
    # 
    # 
    live_variant_caller.write_vcf('output/live_multi_variant_caller.vcf')
    # # live_variant_caller.create_checkpoint('output/checkpoint.pkl')
    #
    # print('total', time.time() - start)


if __name__ == '__main__':
    main()
