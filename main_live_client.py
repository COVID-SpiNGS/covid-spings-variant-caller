import argparse
import socket
import configparser
import logging

logging.basicConfig(filename='vcf_client.log',
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')

parser = argparse.ArgumentParser()
config = configparser.ConfigParser()
config.read('sockets.conf')
HOST = config['BASIC_PARAMS']['HOST']
PORT = int(config['BASIC_PARAMS']['PORT'])
print(PORT)


def _construct_cli():
    #group = parser.add_mutually_exclusive_group()
    #group.add_argument("stop", nargs='?', help="Shutdown server")
    parser.add_argument('processing', help='run or stop', nargs='?', choices=('process', 'write'))
    parser.add_argument('stop', help='run or stop', nargs='?')


def _talk_to_server(action: str, params: str):

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        logging.info(f'Connecting to server under {HOST}:{PORT}...')
        sock.connect((HOST, PORT))
        sock.sendall(b'action')
        sock.close()
        logging.info(f'Closing connection to server under {HOST}:{PORT}...')


if __name__ == '__main__':
    logging.info(f'Welcome... Setting up client')
    _construct_cli()
    args = parser.parse_args()
    action = ''
    params = ''

    print(args)
    if args.stop is not None:
        action = 'stop'

    if args.processing is not None:
        action = 'process'

   # if args.write is not None:
    #    action = 'write'

    logging.info(f'Selected action is {action} with {params}.')

    _talk_to_server(action, params)

