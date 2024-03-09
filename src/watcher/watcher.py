import logging
import os
import sys
import time
from os.path import abspath, dirname

from architecture.client.live_client import VCClient
from config_util import config_io as cio
from config_util.logging import logging as log
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer

log_dir = os.path.join(dirname(dirname(abspath(__file__))), "log")

logging.basicConfig(
    filename=os.path.join(log_dir, "watcher.log"),
    level=logging.DEBUG,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)


class Watcher:
    """
    Class that handles watching a directory for changes
    """

    def __init__(self, directory: str):
        """
        Constructor for Watcher class
        @param directory: relative or absolute path to directory to be watched
        """
        self.directory = directory
        self.recursive = cio.get_watch_recursively()
        self.interval = cio.get_watcher_interval()
        self.client = VCClient(*cio.get_address())
        self.handler = SeqHandler(self.client, cio.get_supported_extensions())
        self.observer = Observer()

    def run(self):
        """
        Function that runs file system watcher in selected directory
        """
        self.observer.schedule(self.handler, self.directory, recursive=self.recursive)
        self.observer.start()
        log.print_and_log(f"Now watching directory {self.directory}", log.INFO)
        try:
            while True:
                time.sleep(self.interval)
        # TODO: Reconsider exception type
        except:
            self.observer.stop()
        self.observer.join()
        log.print_and_log("Watcher terminated.", log.INFO)


class SeqHandler(FileSystemEventHandler):
    """
    Handler class for Watcher, containing logic for handling changes in file system
    """

    def __init__(self, client: VCClient, supported_extensions: list[str]):
        """
        Constructor for SeqHandler
        @param client: instance of VCClient
        @param supported_extensions: supported file extensions to be watched
        """
        self.client = client
        self.supported_extensions = supported_extensions
        self.current_file_size = 0

    def on_any_event(self, event):
        """
        Function to handle any file system event
        @param event: Any file system event
        """
        if not event.is_directory:
            if [
                extension
                for extension in self.supported_extensions
                if event.src_path.endswith(extension)
            ]:
                if event.event_type == "created" or event.event_type == "modified":
                    log.print_and_log(
                        f"Event detector: {event.event_type} in {event.src_path}",
                        log.INFO,
                    )
                    file = event.src_path
                    # stop, process, write
                    self.client.talk_to_server("process", file)
                    # self.client.talk_to_server('write', file)


if __name__ == "__main__":
    path = ""

    if len(sys.argv) > 0:
        path = sys.argv[1]
    else:
        log.print_and_log("No path provided", log.ERROR)

    if os.path.exists(path) and not os.path.isfile(path):
        log.print_and_log(f"Provided path: {path}", log.INFO)
        w = Watcher(path)
        w.run()
    else:
        log.print_and_log(
            f"Provided path {path} does not exist or is a file.", log.ERROR
        )
