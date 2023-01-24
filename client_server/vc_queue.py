import logging
from queue import Queue
from variant_caller.live_variant_caller import LiveVariantCaller
import settings.cio as cio
import configparser
from client_server.vc_exception import VCException

config = configparser.ConfigParser()
config.read('settings.config')

liveVariantCaller = LiveVariantCaller(
    config['VARIANT_CALLER_PARAMS']['REF'],
    cio.get_min_base_quality(),
    cio.get_min_mapping_quality(),
    cio.get_min_total_depth(),
    int(config['VARIANT_CALLER_PARAMS']['minEvidenceDepth']),
    float(config['VARIANT_CALLER_PARAMS']['minEvidenceRatio']),
    int(config['VARIANT_CALLER_PARAMS']['maxVariants'])
)


# try:
# task_queue = VCQueue(queue_size)
# TODO: Reconsider exception type and size
# except VCException:
#   logging.error('Incorrect queue size specified.')

class VCQueue:
    def __init__(self, size: int):
        if 10 >= size >= 1:
            self.size = size
            self.q = Queue(maxsize=self.size)
            logging.info(f'Init queue with size {size}')
        else:
            logging.error(f'Wrong queue size: {size} - Queue size must be in range of (1, 10)!')
            raise VCException(size)

    def put(self, action: str):
        self.q.put(action)

    def put(self, action: (str, str)):
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

    def is_empty(self):
        if self.q.empty():
            return True
        return False

