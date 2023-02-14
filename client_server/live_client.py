import argparse
import socket
import settings.cio as cio
import logging
import os
from pathlib import Path
from os.path import dirname, abspath

log_dir = os.path.join(dirname(dirname(abspath(__file__))), 'log')

logging.basicConfig(filename=os.path.join(log_dir, 'vc_client.log'),
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')

parser = argparse.ArgumentParser()


class VCClient:
    def __init__(self, host, port):

        self.host = host
        self.port = port

    def talk_to_server(self, action: str, params: str):
        payload = bytes(action + ' ' + params, encoding='utf-8')
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
                logging.info(f'Connecting to server under {self.host}:{self.port}...')
                sock.connect((self.host, self.port))
                sock.sendall(payload)
                sock.close()
                logging.info(f'Closing connection to server under {self.host}:{self.port}...')
                print(f'Closing connection to server under {self.host}:{self.port}...')
        except ConnectionRefusedError:
            logging.error(f'Not able to connect to {self.host}:{self.port}. Is server running?')
            print(f'Not able to connect to {self.host}:{self.port}. Is server running?')


def _construct_cli():
    parser.add_argument('--process', help='run or stop', nargs='+')
    parser.add_argument('--write', help='run or stop', nargs='+')
    parser.add_argument('--stop', help='run or stop', nargs='?')


def _params_is_valid(action: str, params: str) -> bool:
    valid = False

    if action.casefold() == 'process':
        if params.endswith('.bam') or params.endswith('.sam') and os.path.isfile(params):
            valid = True

    if action.casefold() == 'write':
        path = Path(params).parent.absolute()
        if params.endswith('.vcf') and os.path.exists(path):
            valid = True

    if action.casefold() == 'stop':
        if params == '':
            valid = True

    return valid


def _run():
    _construct_cli()

    args = parser.parse_args()
    action = ''
    params = ''

    c = VCClient(*cio.get_address())

    if args.stop is not None:
        action = 'stop'

    if args.process is not None:
        action = 'process'
        params = args.process[0]

    if args.write is not None:
        action = 'write'
        params = args.write[0]

    logging.info(f'Selected action is {action} with {params}.')

    if action != '':
        if _params_is_valid(action, params):
            c.talk_to_server(action, params)
        else:
            print(f'{params} is invalid... please make sure path exists.')
            logging.error(f'{params} is invalid... please make sure path exists.')


if __name__ == '__main__':
    logging.info(f'Welcome... Setting up client')
    _run()
