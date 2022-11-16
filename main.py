from variant_caller.live_variant_caller import LiveVariantCaller
from variant_caller.config import minBaseQuality, minMappingQuality, minTotalDepth
import time

def main():
    start = time.time()

    minEvidenceDepth = 5
    minEvidenceRatio = 0.0
    maxVariants = 1


    liveVariantCaller = LiveVariantCaller(
        'input/reference.fasta',
        minBaseQuality,
        minMappingQuality,
        minTotalDepth,
        minEvidenceDepth,
        minEvidenceRatio,
        maxVariants
    )

    # liveVariantCaller.process_bam('input/input.bam')


    start_sub = time.time()
    liveVariantCaller.process_bam('input/input-1/input.bam')
    print('1', time.time() - start_sub)
  

    start_sub = time.time()
    liveVariantCaller.process_bam('input/input-2/input.bam')
    print('2', time.time() - start_sub)


    start_sub = time.time()
    liveVariantCaller.process_bam('input/input-3/input.bam')
    print('3', time.time() - start_sub)

    start_sub = time.time()
    liveVariantCaller.process_bam('input/input-4/input.bam')
    print('4', time.time() - start_sub)
    
    start_sub = time.time()
    liveVariantCaller.process_bam('input/input-5/input.bam')
    print('5', time.time() - start_sub)
    
    start_sub = time.time()
    liveVariantCaller.process_bam('input/input-6/input.bam')
    print('6', time.time() - start_sub)
    
    start_sub = time.time()
    liveVariantCaller.process_bam('input/input-7/input.bam')
    print('7', time.time() - start_sub)
    
    start_sub = time.time()
    liveVariantCaller.process_bam('input/input-8/input.bam')
    print('8', time.time() - start_sub)
    
    start_sub = time.time()
    liveVariantCaller.process_bam('input/input-9/input.bam')
    print('9', time.time() - start_sub)
    
    start_sub = time.time()
    liveVariantCaller.process_bam('input/input-10/input.bam')
    print('10', time.time() - start_sub)


    liveVariantCaller.write_vcf('output/live_multi_variant_caller.vcf')
    # liveVariantCaller.create_checkpoint('output/checkpoint.pkl')

    print('total', time.time() - start)
  

if __name__=='__main__':
    main()