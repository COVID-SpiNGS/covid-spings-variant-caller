import logging
from queue import Queue
from variant_caller.live_variant_caller import LiveVariantCaller
from variant_caller.config import minBaseQuality, minMappingQuality, minTotalDepth
import configparser
from vc_exception import VCException

config = configparser.ConfigParser()
config.read('settings.config')

liveVariantCaller = LiveVariantCaller(
    config['VARIANT_CALLER_PARAMS']['REF'],
    minBaseQuality,
    minMappingQuality,
    minTotalDepth,
    int(config['VARIANT_CALLER_PARAMS']['minEvidenceDepth']),
    float(config['VARIANT_CALLER_PARAMS']['minEvidenceRatio']),
    int(config['VARIANT_CALLER_PARAMS']['maxVariants'])
)


class VCQueue:
    q = Queue(maxsize=1)

    def __init__(self, size: int):
        if 10 >= size >= 1:
            self.size = size
            self.q = Queue(maxsize=self.size)
            logging.info(f'Init queue with size {size}')
        else:
            raise VCException(size)

    def put(self, action: str):
        self.q.put(action)

    def process(self):
        if not self.q.empty():
            (action, path) = self.q.get()

            if action == 'process':
                self._process_bam(path)

            elif action == 'write':
                self._write_vcf(path)

    def _write_vcf(self, path: str):
        logging.info(f'Writing VCF to {path}')
        liveVariantCaller.write_vcf(path)

    def _process_bam(self, path: str):
        logging.info(f'Processing BAM with path {path}')
        liveVariantCaller.process_bam(path)
