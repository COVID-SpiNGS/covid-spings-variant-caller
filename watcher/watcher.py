from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from client_server.live_client import VCClient
import settings.cio as cio
import logging
import time
import sys
import os

logging.basicConfig(filename='../log/watcher.log',
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')


class Watcher:

    def __init__(self, directory):
        self.directory = directory
        self.recursive = cio.get_watch_recursively()
        self.interval = cio.get_watcher_interval()
        self.client = VCClient(*cio.get_address())
        self.handler = SeqHandler(cio.get_supported_extensions())
        self.observer = Observer()


    def run(self):
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

    def __init__(self, supported_extensions):
        self.supported_extensions = supported_extensions

    def on_any_event(self, event):
        file_extension = os.path.splitext(event.src_path)
        if not event.is_directory and file_extension[-1] == '.txt':
            if event.event_type == 'created' or event.event_type == 'modified':
                logging.info(f'Event detector: {event.event_type} in {event.src_path}')
                print(f'Event detector: {event.event_type} in {event.src_path}')

                if path.endswith(self.supported_extensions):
                    print('lol')

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
