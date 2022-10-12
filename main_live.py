import sys
import socket
import time
import threading
import daemon
import logging
import argparse

logging.basicConfig(filename='vcf.log',
level=logging.DEBUG,
format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')


HOST = '127.0.0.1'  # Standard loopback interface address (localhost)
PORT = 65432  # Port to listen on (non-privileged ports are > 1023)

parser = argparse.ArgumentParser()


def handle_client(sock):
    with sock.makefile() as f:
        sock.close()
        for line in f:
            f.writeline(line)


def serve_forever():
    server = socket.socket()
    server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    server.bind((HOST, PORT))
    server.listen(1)
    while True:
        conn, address = server.accept()
        thread = threading.Thread(target=handle_client, args=[conn])
        thread.daemon = True
        thread.start()



def construct_cli():
    group = parser.add_mutually_exclusive_group()
    group.add_argument("start", nargs='?', help="echo the string you use here")
    group.add_argument("stop", nargs='?', help="echo the string you use here")
    group.add_argument("process", nargs='?', help="echo the string you use here")
    group.add_argument("write", nargs='?', help="echo the string you use here")


def _process_bam(path: str):
    logging.info('Processing BAM...', path)
    return True

def _write_vcf(path: str):
    return True

def _shutdown_gracefully(sock):
    logging.info('Stopping server in 10 seconds...')
    time.sleep(10)
    sock.shutdown(socket.SHUT_RDWR)
    sock.close()
    return True

if __name__ == '__main__':

    sock_running = False
    construct_cli()

    args = parser.parse_args()

    if args.start is not None and sock_running is False:
        #sock_running = True
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.bind((HOST, PORT))
            sock.listen()
            logging.info(f'Running now under {HOST}:{PORT}...')

            while True:
                connection, address = sock.accept()
                with connection:
                    data = connection.recv(1024)
                    logging.info(f"Received {data!r}")
                    action = b''
                    params = b''

                    if b' ' in data:
                        (action, params) = data.split(b' ')
                    else:
                        action = data
                        params = b''

                    if action == b'stop':
                        _shutdown_gracefully(sock)
                    elif action == b'process':

                        _process_bam(params)
                    elif action == b'write':
                        logging.info('Write VCF...', params)
                        _write_vcf(params)
                    else:
                        logging.info("NO SUCH ACTION")
                        # TODO: Throw Illegal Action Exception maybe ?
    #with daemon.DaemonContext():
    #    logging.info("LOL")
    #serve_forever()
