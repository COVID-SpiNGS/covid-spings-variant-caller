import logging
import os
from queue import Queue
from variant_caller.live_variant_caller import LiveVariantCaller
import settings.cio as cio
from client_server.vc_exception import VCException
from os.path import dirname, abspath


# try:
# task_queue = VCQueue(queue_size)
# TODO: Reconsider exception type and size
# except VCException:
#   logging.error('Incorrect queue size specified.')

log_dir = os.path.join(dirname(dirname(abspath(__file__))), 'log')

logging.basicConfig(filename=os.path.join(log_dir, 'vc_server.log'),
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s - VC_Q |')

class VCQueue:
    def __init__(self, size: int):
        if 10 >= size >= 1:
            self.size = size
            self.q = Queue(maxsize=self.size)
            self.live_variant_caller = LiveVariantCaller(
                '../' + cio.get_reference(),
                cio.get_min_base_quality(),
                cio.get_min_mapping_quality(),
                cio.get_min_total_depth(),
                cio.get_min_evidence_depth(),
                cio.get_min_evidence_ratio(),
                cio.get_max_variants()
            )
            self.temp_dir = cio.get_temp_dir()
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
        self.live_variant_caller.write_vcf(path)

    def _process_bam(self, path: str):
        logging.info(f'Processing BAM with path {path}')
        checkpoint = self.temp_dir + '/' + os.path.basename(path)
        if os.path.exists(checkpoint):
            self.live_variant_caller.load_checkpoint(checkpoint)

        self.live_variant_caller.process_bam(path)
        self.live_variant_caller.create_checkpoint(checkpoint)

    def is_empty(self):
        return self.q.empty()
