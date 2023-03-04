import logging
from time import strftime, localtime

DEBUG = 'debug'
ERROR = 'error'
INFO = 'info'
WARNING = 'warning'


def print_and_log(text: str, log_type: str):
    timestamp = strftime('[%Y-%m-%d %H:%M:%S]', localtime())
    print(f'{timestamp} {text}')

    if log_type == 'debug':
        logging.debug(text)
    if log_type == 'error':
        logging.error(text)
    if log_type == 'info':
        logging.info(text)
    if log_type == 'warning':
        logging.warning(text)
