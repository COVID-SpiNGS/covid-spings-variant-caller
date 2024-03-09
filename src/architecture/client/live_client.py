import argparse
import socket
import config_util.cio as cio
import logging
import os
from pathlib import Path
from os.path import dirname, abspath
import config_util.logging as log

log_dir = os.path.join(dirname(dirname(abspath(__file__))), "log")

logging.basicConfig(
    filename=os.path.join(log_dir, "vc_client.log"),
    level=logging.DEBUG,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)

parser = argparse.ArgumentParser()


class VCClient:
    """
    Class of Variant Caller client

    This class represents a client that connects to a Variant Caller service. It is used to interact with the VC server.

    Attributes:
        host (str): The host address of the Variant Caller server.
        port (int): The port number on which the Variant Caller server is listening.

    """

    def __init__(self, host: str, port: int):
        """
        Constructor for client
        @param host: host address
        @param port: host port
        """
        self.host = host
        self.port = port

    def talk_to_server(self, action: str, path: str):
        """
        Function to send corresponding message to server
        @param action: action to be executed
        @param path: path to file for action
        """
        payload = bytes(action + " " + path, encoding="utf-8")
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
                logging.info(f"Connecting to server under {self.host}:{self.port}...")
                sock.connect((self.host, self.port))
                sock.sendall(payload)
                sock.close()
                log.print_and_log(
                    f"Closing connection to server under {self.host}:{self.port}...",
                    log.INFO,
                )

        except ConnectionRefusedError:
            log.print_and_log(
                f"Not able to connect to {self.host}:{self.port}. Is server running?",
                log.ERROR,
            )


def _construct_cli():
    """
    Function that creates command line interface
    """
    parser.add_argument("-p", "--process", help="run or stop", nargs="+")
    parser.add_argument("-w", "--write", help="run or stop", nargs="+")
    parser.add_argument("-st", "--stop", help="run or stop", action="store_true")


def _params_is_valid(action: str, path: str) -> bool:
    """
    Function that validates input form CLI
    @param action: action to be executedd
    @param path: path provided by CLI
    @return: bool whether params are valid (path exists,...)
    """
    valid = False

    if action.casefold() == "process":
        if path.endswith(".bam") and os.path.isfile(path):
            valid = True

    if action.casefold() == "write":
        path = Path(path).parent.absolute()
        if path.endswith(".vcf") and os.path.exists(path):
            valid = True

    if action.casefold() == "stop":
        if path == "":
            valid = True

    return valid


def _run():
    """
    Function that runs client with params from CLI
    """
    _construct_cli()

    args = parser.parse_args()
    action = ""
    path = ""

    c = VCClient(*cio.get_address())

    if args.stop is not None:
        action = "stop"

    if args.process is not None:
        action = "process"
        path = args.process[0]

    if args.write is not None:
        action = "write"
        path = args.write[0]

    logging.info(f"Selected action is {action} with {path}.")

    if action != "":
        if _params_is_valid(action, path):
            c.talk_to_server(action, path)
        else:
            log.print_and_log(
                f"{path} is invalid... please make sure path exists.", log.ERROR
            )


if __name__ == "__main__":
    logging.info("Welcome... Setting up client")
    _run()
