import multiprocessing
import socket
import time
import threading
import daemon
from client_server.vc_exception import VCException
from client_server.vc_queue import VCQueue
import logging
import configparser

logging.basicConfig(filename='../log/vc_server.log',
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')


class VCServer:

    def __init__(self, host, port, queue_size):
        self.host = host
        self.port = port
        self.queue_size = queue_size
        self.task_queue = VCQueue(self.queue_size)

    def run(self):

        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.bind((self.host, self.port))
            sock.listen()
            logging.info(f'Running now under {self.host}:{self.port}...')
            print(f'Running now under {self.host}:{self.port}...')

            while True:
                connection, address = sock.accept()
                with connection:
                    data = connection.recv(1024)
                    logging.info(f"Received {data!r}")
                    recv_data = data.decode('utf-8').split(' ')

                    if recv_data[0] == 'stop':
                        self._shutdown_gracefully(sock)
                    elif recv_data[0] == 'process' or recv_data[0] == 'write':

                        self.task_queue.put((recv_data[0], recv_data[1]))

                    else:
                        logging.error(f'No such action: {recv_data[0]}')

                    while not self.task_queue.is_empty():
                        self.task_queue.process()

    def _shutdown_gracefully(self, sock):
        logging.info('Stopping server in 10 seconds...')
        print('Stopping server in 10 seconds...')
        time.sleep(10)
        sock.shutdown(socket.SHUT_RDWR)
        sock.close()
        return True


if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read('settings.config')
    host = config['BASIC_PARAMS']['HOST']
    port = int(config['BASIC_PARAMS']['PORT'])
    queue_size = int(config['BASIC_PARAMS']['QUEUE_SIZE'])
    server = VCServer(host, port, queue_size)
    server.run()
