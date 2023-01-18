from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import client_server
import logging
import time
import argparse
import sys
import os

logging.basicConfig(filename='../log/watcher.log',
                    level=logging.DEBUG,
                    format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')


class Watcher:

    def __init__(self, directory, handler=FileSystemEventHandler()):
        self.observer = Observer()
        self.handler = handler
        self.directory = directory

    def run(self):
        self.observer.schedule(
            self.handler, self.directory, recursive=True)
        self.observer.start()
        logging.info(f'Watcher running in {self.directory}')
        print(f'Now watching directory {self.directory}')
        try:
            while True:
                time.sleep(1)
        # TODO: Reconsider exception type
        except:
            self.observer.stop()
        self.observer.join()
        print('Watcher terminated.')
        logging.info('Watcher terminated.')


class SeqHandler(FileSystemEventHandler):

    def on_any_event(self, event):
        file_extension = os.path.splitext(event.src_path)
        if not event.is_directory and file_extension[-1] == '.txt':
            if event.event_type == 'created' or event.event_type == 'modified':
                logging.info(f'Event detector: {event.event_type} in {event.src_path}')
                print(f'Event detector: {event.event_type} in {event.src_path}')  # Your code here


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
        w = Watcher(path, SeqHandler())
        w.run()
    else:
        logging.error(f'Path {path} does not exist or is a file.')
        print(f'The provided path - {path} - does not exist or is a file. Path must be a directory!')
