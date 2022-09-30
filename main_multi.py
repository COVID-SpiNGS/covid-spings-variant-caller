from variant_caller.live_multi_variant_caller import LiveMultiVariantCaller
from config import minBaseQuality, minMappingQuality, minTotalDepth

def main():
    minEvidenceDepth = 5
    minEvidenceRatio = 0.10


    liveVariantCaller = LiveMultiVariantCaller(
        'input/reference.fasta',
        minBaseQuality,
        minMappingQuality,
        minTotalDepth,
        minEvidenceDepth,
        minEvidenceRatio
    )


    liveVariantCaller.process_bam('input/input-1/input.bam')
    liveVariantCaller.process_bam('input/input-2/input.bam')
    liveVariantCaller.process_bam('input/input-3/input.bam')
    liveVariantCaller.process_bam('input/input-4/input.bam')
    liveVariantCaller.process_bam('input/input-5/input.bam')
    liveVariantCaller.process_bam('input/input-6/input.bam')
    liveVariantCaller.process_bam('input/input-7/input.bam')
    liveVariantCaller.process_bam('input/input-8/input.bam')
    liveVariantCaller.process_bam('input/input-9/input.bam')
    liveVariantCaller.process_bam('input/input-10/input.bam')
    
    liveVariantCaller.write_vcf('output/live_multi_variant_caller.vcf')
    liveVariantCaller.create_checkpoint('output/checkpoint.pkl')
  

if __name__=='__main__':
    main()