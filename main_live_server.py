import multiprocessing
import socket
import time
import threading
import daemon
from vc_exception import VCException
from vc_queue import VCQueue
import logging
import configparser

logging.basicConfig(filename='vc_server.log',
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')

config = configparser.ConfigParser()
config.read('settings.config')
HOST = config['BASIC_PARAMS']['HOST']
PORT = int(config['BASIC_PARAMS']['PORT'])
queue_size = int(config['BASIC_PARAMS']['QUEUE_SIZE'])
task_queue = VCQueue(queue_size)


def _shutdown_gracefully(sock):
    logging.info('Stopping server in 10 seconds...')
    time.sleep(10)
    sock.shutdown(socket.SHUT_RDWR)
    sock.close()
    return True


def _run():
    task_queue = VCQueue(queue_size)

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind((HOST, PORT))
        sock.listen()
        logging.info(f'Running now under {HOST}:{PORT}...')

        while True:
            connection, address = sock.accept()
            with connection:
                data = connection.recv(1024)
                logging.info(f"Received {data!r}")
                recv_data = data.decode('utf-8').split(' ')

                if recv_data[0] == 'stop':
                    _shutdown_gracefully(sock)
                elif recv_data[0] == 'process' or recv_data[0] == 'write':

                    task_queue.put((recv_data[0], recv_data[1]))

                else:
                    logging.error(f'No such action: {recv_data[0]}')

                while not task_queue.is_empty():
                    task_queue.process()


# with daemon.DaemonContext():
#    logging.info("LOL")
# serve_forever()


if __name__ == '__main__':
    _run()
