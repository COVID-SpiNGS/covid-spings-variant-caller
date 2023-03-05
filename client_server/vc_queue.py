import logging
import os
import threading
from queue import Queue
from variant_caller.live_variant_caller import LiveVariantCaller
import config_util.cio as cio
from client_server.vc_exception import VCException
from os.path import dirname, abspath
import pysam as pys
import config_util.logging as log

# try:
# task_queue = VCQueue(queue_size)
# TODO: Reconsider exception type and size
# except VCException:
#   logging.error('Incorrect queue size specified.')

log_dir = os.path.join(dirname(dirname(abspath(__file__))), 'log')
logging.basicConfig(filename=os.path.join(log_dir, 'vc_server.log'),
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')


def _run_samtools(sam_file_path: str, abs_path: str, bam_file_name: str):
    """
    Function that runs samtools sort and index when .sam file changes are monitored
    @param sam_file_path: path to sam file path to be processed
    @param abs_path: absolute path to .sam file folder
    @param: bam_file_name: .bam file name to be used for output (must be equal to .sam file as per convention)
    """
    # Based on command: (['samtools', 'sort', '-O', 'bam', '-o', 'sorted_test_seq.bam', 'test_seq.sam'])

    try:
        pys.sort('-O', 'bam', '-o', os.path.join(abs_path, bam_file_name), sam_file_path)
        pys.index(os.path.join(abs_path, bam_file_name), os.path.join(abs_path, bam_file_name + '.bai'))
    except Exception as e:
        log.print_and_log(e, log.ERROR)
        return


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
            self.output_dir = cio.get_output_dir()
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
        log.print_and_log(f'Added {action} to queue', log.INFO)
        log.print_and_log(f'Current queue size: {self.q.qsize()}', log.INFO)

    def process(self):
        """
        Implementation of behaviour when element is retrieved from queue
        """

        if not self.q.empty():
            (action, path) = self.q.get()
            log.print_and_log(f'Queue size atm is {self.q.qsize()}', log.DEBUG)
            log.print_and_log(f'Current action is: {action}', log.DEBUG)

            daemon = threading.Thread(name='daemon_vc', target=self._process_bam, args=(path, ))

            #if action == 'process':
                #self._process_bam(path)

            if action == 'write':
                daemon = threading.Thread(name='daemon_vc', target=self._write_vcf, args=(path, ))
                self._write_vcf(path)

            self.current_size -= 1

            daemon.daemon = True
            daemon.start()
            #print(f'Queue size atm is {self.q.qsize()} - {self.current_size}')

    def _write_vcf(self, path: str):
        """
        Function acting as wrapper for variant caller's function to write VCF report
        @param path: path to vcf file
        """
        vcf_path = path.split('.bam')[0] + '.vcf'
        log.print_and_log(f'Writing VCF to {vcf_path}', log.INFO)
        self.live_variant_caller.write_vcf(vcf_path)

    def _process_bam(self, file_path: str):
        """
        Function acting as wrapper for variant caller's function to process BAM file
        @param path: path to BAM file
        """
        logging.info(f'Processing BAM with path {file_path}')
        abs_path = os.path.dirname(file_path)
        file_name = os.path.basename(file_path)
        if file_name.endswith(cio.SAM):
            file_name = file_name.split(cio.SAM)[0] + cio.BAM
        _run_samtools(file_path, abs_path, file_name)
        checkpoint = os.path.join(self.temp_dir, file_name + cio.get_temp_file_extension())
        index_file = file_name + cio.BAI

        if os.path.exists(file_path) and os.path.exists(os.path.join(abs_path, index_file)):
            if os.path.exists(checkpoint):
                log.print_and_log(f'Checkpoint for {file_name} found', log.DEBUG)
                self.live_variant_caller.load_checkpoint(checkpoint)

            self.live_variant_caller.process_bam(os.path.join(abs_path, file_name))
            self.live_variant_caller.create_checkpoint(checkpoint)
            self.live_variant_caller.write_vcf(os.path.join(self.output_dir, file_name + cio.VCF))
        else:
            log.print_and_log(f'{file_name} or {index_file} or do not exist', log.ERROR)

    def length(self) -> int:
        """
        Wrapper function for present number of elements in queue
        @return: number of current elements in queue
        """
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
