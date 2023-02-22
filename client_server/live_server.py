import multiprocessing
import socket
import time
import threading
import daemon
import os
import config.cio as cio
from client_server.vc_exception import VCException
from os.path import dirname, abspath
from client_server.vc_queue import VCQueue
import logging

log_dir = os.path.join(dirname(dirname(abspath(__file__))), 'log')

logging.basicConfig(filename=os.path.join(log_dir, 'vc_server.log'),
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')


class VCServer:
    """

    """

    def __init__(self):
        """

        """
        host, port = cio.get_address()
        self.host = host
        self.port = port
        self.queue_size = cio.get_queue_size()
        self.task_queue = VCQueue(self.queue_size)

    def run(self):
        """

        """
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.bind((self.host, self.port))
            sock.listen()
            logging.info(f'Running now under {self.host}:{self.port}...')
            print(f'Running now under {self.host}:{self.port}...')

            while True:
                connection, address = sock.accept()
                with connection:
                    data = connection.recv(1024)
                    logging.info(f'Received {data!r}')
                    print(f'Received {data!r}')

                    recv_data = data.decode('utf-8').split(' ')

                    if recv_data[0] == 'stop':
                        ret = self._shutdown_gracefully(sock)
                    elif recv_data[0] == 'process' or recv_data[0] == 'write':
                        logging.info(f'Received {data[0]} with argument {data[1]}')
                        if self.task_queue.length() < self.queue_size:
                            self.task_queue.put((recv_data[0], recv_data[1]))
                        else:
                            pass
                            #TODO: ADD MAX QUEUE ERROR

                    else:
                        logging.error(f'No such action: {recv_data[0]}')
                        print(f'No such action: {recv_data[0]}')

                    while not self.task_queue.is_empty():
                        self.task_queue.process()
                        self.task_queue.join()


    def _shutdown_gracefully(self, sock):
        """

        @param sock:
        @return:
        """
        try:
            logging.info('Stopping server in 10 seconds...')
            print('Stopping server in 10 seconds...')
            time.sleep(10)
            sock.shutdown(socket.SHUT_RDWR)
            sock.close()
            return 0
        except:
            return -1


if __name__ == '__main__':
    server = VCServer()
    server.run()
