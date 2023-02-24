import logging
import os
import time
from queue import Queue
from variant_caller.live_variant_caller import LiveVariantCaller
import config.cio as cio
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
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')


class VCQueue:
    """
    Wrapper class around multiprocessing.queue.Queue with additional business logic
    """
    def __init__(self, size: int):
        """
        Constructor for VCQueue
        @param size: queue size (may not exceed 10)
        """
        if cio.get_max_queue_size() >= size >= cio.get_min_queue_size():
            self.size = size
            self.q = Queue(maxsize=self.size)
            self.current_size = 0
            self.live_variant_caller = LiveVariantCaller(
                os.path.join(dirname(dirname(abspath(__file__))), cio.get_reference()),
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
        """
        Add element to queue
        @param action: action to be performed as string
        """
        self.q.put(action)
        self.current_size += 1


    def put(self, action: (str, str)):
        """
        Add action  to queue
        @param action: tuple (action, file path) to be added to queue
        """
        self.q.put(action)
        self.current_size += 1
        print(f'Added {action} to queue')
        #print(f' 2 Queue size atm is {self.q.qsize()}')

    def process(self):
        """

        """
        if not self.q.empty():
            (action, path) = self.q.get()
            logging.debug(f'Queue size atm is {self.q.qsize()}')
            #print(f'Queue size atm is {self.q.qsize()}')
            #time.sleep(10)
            print(f'Current action is: {action}')
            if action == 'process':
                self._process_bam(path)

            elif action == 'write':
                self._write_vcf(path)

            self.current_size -= 1
            self.q.task_done()
            print(f'Queue size atm is {self.q.qsize()} - {self.current_size}')

    def _write_vcf(self, path: str):
        """
        @param path:
        """
        vcf_path = path.split('.bam')[0] +'.vcf'
        logging.info(f'Writing VCF to {vcf_path}')
        print(f'Writing VCF to {vcf_path}')
        self.live_variant_caller.write_vcf(vcf_path)

    def _process_bam(self, path: str):
        """

        @param path:
        """
        logging.info(f'Processing BAM with path {path}')
        basename = os.path.basename(path)
        checkpoint = os.path.join(self.temp_dir, basename + '.pkl')
        index_file = path + '.bai'

        if os.path.exists(path) and os.path.exists(index_file):
            if os.path.exists(checkpoint):
                print(f'Checkpoint for {basename} found')
                logging.debug(f'Checkpoint for {basename} found')
                self.live_variant_caller.load_checkpoint(checkpoint)

            print(f'running under: {path}')
            self.live_variant_caller.process_bam(path)
            self.live_variant_caller.create_checkpoint(checkpoint)
        else:
            logging.error(f'{path} or {index_file} or do not exist')
            print(f'{path} or {index_file} do not exist')

    def length(self) -> int:
        """
        Wrapper function for present number of elements in queue
        @return: number of current elements in queue
        """
        logging.info(f'Queue size - {self.q.qsize()}')
        return self.q.qsize()

    def is_empty(self) -> bool:
        """
        Wrapper function to check whether queue is empty
        @return: bool whether queue is empty
        """
        return self.q.empty()

    def join(self):
        """
        Wrapper for queue.join() - wait until task is done
        """
        self.q.join()
