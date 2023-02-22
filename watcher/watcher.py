from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from client_server.live_client import VCClient
import config.cio as cio
import logging
import time
import sys
import os
from os.path import dirname, abspath

log_dir = os.path.join(dirname(dirname(abspath(__file__))), 'log')

logging.basicConfig(filename=os.path.join(log_dir, 'watcher.log'),
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')


class Watcher:
    """

    """

    def __init__(self, directory):
        """

        @param directory:
        """
        self.directory = directory
        self.recursive = cio.get_watch_recursively()
        self.interval = cio.get_watcher_interval()
        self.client = VCClient(*cio.get_address())
        self.handler = SeqHandler(self.client, cio.get_supported_extensions())
        self.observer = Observer()

    def run(self):
        """

        @return:
        """
        self.observer.schedule(
            self.handler, self.directory, recursive=self.recursive)
        self.observer.start()
        logging.info(f'Watcher running in {self.directory}')
        print(f'Now watching directory {self.directory}')
        try:
            while True:
                time.sleep(self.interval)
        # TODO: Reconsider exception type
        except:
            self.observer.stop()
        self.observer.join()
        print('Watcher terminated.')
        logging.info('Watcher terminated.')


class SeqHandler(FileSystemEventHandler):
    """

    """

    def __init__(self, client, supported_extensions):
        """

        @param client:
        @param supported_extensions:
        """
        self.client = client
        self.supported_extensions = supported_extensions

    def on_any_event(self, event):
        """

        @param event:
        @return:
        """
        if not event.is_directory:
            if [extension for extension in self.supported_extensions if event.src_path.endswith(extension)]:
                if event.event_type == 'created' or event.event_type == 'modified':
                    logging.info(f'Event detector: {event.event_type} in {event.src_path}')
                    print(f'Event detector: {event.event_type} in {event.src_path}')
                    file = event.src_path
                    # stop, process, write
                    self.client.talk_to_server('process', file)
                    self.client.talk_to_server('write', file)

    def on_modified(self, event):
        pass

    def on_deleted(self, event):
        pass

    def on_created(self, event):
        pass

    def on_moved(self, event):
        pass


if __name__ == '__main__':

    path = ''

    if len(sys.argv) > 0:
        path = sys.argv[1]
    else:
        logging.error(f'No path provided')
        print(f'Please provide path to be watched!')

    if os.path.exists(path) and not os.path.isfile(path):
        logging.info(f'Provided path: {path}')
        print(f'Provided path: {path}')
        w = Watcher(path)
        w.run()
    else:
        logging.error(f'Path {path} does not exist or is a file.')
        print(f'The provided path - {path} - does not exist or is a file. Path must be a directory!')
