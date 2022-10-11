import sys
import socket
import time
import threading
import daemon
import argparse


HOST = '127.0.0.1'  # Standard loopback interface address (localhost)
PORT = 65432  # Port to listen on (non-privileged ports are > 1023)


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



parser = argparse.ArgumentParser()

def construct_cli():
    parser.add_argument("start", help="echo the string you use here")
    parser.add_argument("stop", help="echo the string you use here")
    parser.add_argument("process", help="echo the string you use here")
    parser.add_argument("write", help="echo the string you use here")


if __name__ == '__main__':

    construct_cli()

    args = parser.parse_args()
    print(args.echo)

    #with daemon.DaemonContext():
    #    print("LOL")
    #serve_forever()
