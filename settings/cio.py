import configparser

config = configparser.ConfigParser()
config.read('settings.config')


# Basic params
def get_address() -> (str, int):
    return '127.0.0.1', 65342#config['BASIC_PARAMS']['HOST'], int(config['BASIC_PARAMS']['PORT'])


def get_queue_size() -> int:
    return int(config['BASIC_PARAMS']['QUEUE_SIZE'])


# Variant Caller Params


# Watcher Params

def get_watcher_interval() -> int:
    return int(config['WATCHER_PARAMS']['WATCHER_INTERVAL'])


def get_watch_recursively() -> bool:
    return bool(config['WATCHER_PARAMS']['WATCH_RECURSIVELY'])


def get_supported_extensions() -> str:
    return config['WATCHER_PARAMS']['SUPPORTED_EXTENSIONS']
