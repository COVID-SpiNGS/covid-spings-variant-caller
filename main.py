from variant_caller.live_variant_caller import LiveVariantCaller
from config import minBaseQuality, minMappingQuality, minCandidatesDepth, minTotalDepth

def main():
    liveVariantCaller = LiveVariantCaller(
        'input/reference.fasta',
        minBaseQuality,
        minMappingQuality,
        minTotalDepth,
        minCandidatesDepth
    )

    liveVariantCaller.process_bam('input/input.bam')
    liveVariantCaller.process_bam('input/input.bam')
    liveVariantCaller.write_vcf('output/live_variant_caller.vcf')
  

if __name__=='__main__':
    main()