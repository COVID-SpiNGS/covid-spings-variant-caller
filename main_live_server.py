import os
import socket
import time
import threading
import daemon
import logging
import configparser

logging.basicConfig(filename='vcf_server.log',
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')

config = configparser.ConfigParser()
config.read('sockets.conf')
HOST = config['BASIC_PARAMS']['HOST']
PORT = int(config['BASIC_PARAMS']['PORT'])


def _process_bam(path: str):
    logging.info('Processing BAM...', path)
    return True


def _write_vcf(path: str):
    logging.info('Write VCF...', path)
    return True


def _path_is_valid(action: str, path: str) -> bool:
    valid = False

    if action.casefold() == 'process':
        if path.endswith(b'.bam') and os.path.isfile(path):
            valid = True

    if action.casefold() == 'write':
        if path.endswith(b'.vcf') and os.path.isfile(path):
            valid = True

    return valid


def _shutdown_gracefully(sock):
    logging.info('Stopping server in 10 seconds...')
    time.sleep(10)
    sock.shutdown(socket.SHUT_RDWR)
    sock.close()
    return True


def _run():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind((HOST, PORT))
        sock.listen()
        logging.info(f'Running now under {HOST}:{PORT}...')

        while True:
            connection, address = sock.accept()
            with connection:
                data = connection.recv(1024)
                logging.info(f"Received {data!r}")
                (action, params) = (b'', b'')

                if b' ' in data:
                    (action, params) = data.split(b' ')
                else:
                    action = data
                    params = b''

                if action == b'stop':
                    _shutdown_gracefully(sock)
                elif action == b'process':
                    if _path_is_valid(action, params):
                        _process_bam(params)
                elif action == b'write':
                    if _path_is_valid(action, params):
                        _write_vcf(params)
                else:
                    logging.info("NO SUCH ACTION")
                    # TODO: Throw Illegal Action Exception maybe ?


# with daemon.DaemonContext():
#    logging.info("LOL")
# serve_forever()


if __name__ == '__main__':
    _run()
