import logging
import os
import socket
import time
from os.path import abspath, dirname


from src.config_util import logging as log
from src.architecture.vc_queue import VCQueue
from src.config_util import config_io as cio

log_dir = os.path.join(dirname(dirname(abspath(__file__))), "log")

logging.basicConfig(
    filename=os.path.join(log_dir, "vc_server.log"),
    level=logging.DEBUG,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)


class VCServer:
    """
    Class of Variant Caller server

    Attributes:
            host (str): The host address obtained from 'config_io'.
            port (int): The port number obtained from 'config_io'.
            queue_size (int): The size of the task queue, obtained from 'config_io'.
            task_queue (VCQueue): An instance of 'VCQueue' initialized with the 'queue_size'.

    """

    def __init__(self):
        """
        Constructor for VC server
        """
        host, port = cio.get_address()
        self.host = host
        self.port = port
        self.queue_size = cio.get_queue_size()
        self.task_queue = VCQueue(self.queue_size)

    def run(self):
        """
        Function that runs server on given host and port
        """
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.bind((self.host, self.port))
            sock.listen()
            log.print_and_log(f"Running now under {self.host}:{self.port}...", log.INFO)

            while True:
                connection, address = sock.accept()
                with connection:
                    data = connection.recv(1024)
                    log.print_and_log(f"Received {data!r}", log.INFO)

                    recv_data = data.decode("utf-8").split(" ")

                    if recv_data[0] == "stop":
                        ret = self._shutdown_gracefully(sock)
                        break

                    elif recv_data[0] == "process" or recv_data[0] == "write":
                        logging.info(f"Received {data[0]} with argument {data[1]}")
                        if self.task_queue.length() < self.queue_size:
                            self.task_queue.put((recv_data[0], recv_data[1]))
                        else:
                            pass
                            # TODO: ADD MAX QUEUE ERROR

                    else:
                        log.print_and_log(f"No such action: {recv_data[0]}", log.ERROR)

                    while not self.task_queue.is_empty():
                        self.task_queue.process()
                        # self.task_queue.join()

    def _shutdown_gracefully(self, sock):
        """
        Function that shuts down server socket upon corresponding message
        @param sock: socket to be shutdown
        """
        try:
            log.print_and_log("Stopping server in 10 seconds...", log.INFO)
            time.sleep(10)
            sock.shutdown(socket.SHUT_RDWR)
            sock.close()
            return 0
        except:
            return -1


if __name__ == "__main__":
    server = VCServer()
    server.run()
