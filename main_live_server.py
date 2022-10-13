import socket
import time
import threading
import daemon
import logging
import configparser
from variant_caller.live_variant_caller import LiveVariantCaller
from variant_caller.config import minBaseQuality, minMappingQuality, minTotalDepth

logging.basicConfig(filename='vcf_server.log',
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')

config = configparser.ConfigParser()
config.read('settings.config')
HOST = config['BASIC_PARAMS']['HOST']
PORT = int(config['BASIC_PARAMS']['PORT'])

liveVariantCaller = LiveVariantCaller(
    config['VARIANT_CALLER_PARAMS']['REF'],
    minBaseQuality,
    minMappingQuality,
    minTotalDepth,
    int(config['VARIANT_CALLER_PARAMS']['minEvidenceDepth']),
    float(config['VARIANT_CALLER_PARAMS']['minEvidenceRatio']),
    int(config['VARIANT_CALLER_PARAMS']['maxVariants'])
)


def _process_bam(path: str):
    logging.info(f'Processing BAM with path {path}')
    liveVariantCaller.process_bam(path)


def _write_vcf(path: str):
    logging.info(f'Writing VCF to {path}')
    liveVariantCaller.write_vcf(path)


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
                recv_data = data.decode('utf-8').split(' ')

                if recv_data[0] == 'stop':
                    _shutdown_gracefully(sock)
                elif recv_data[0] == 'process':
                    _process_bam(recv_data[1])
                elif recv_data[0] == 'write':
                    _write_vcf(recv_data[1])
                else:
                    logging.error(f'No such action: {recv_data[0]}')


# with daemon.DaemonContext():
#    logging.info("LOL")
# serve_forever()


if __name__ == '__main__':
    _run()
