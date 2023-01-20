import configparser


def get_configparser():
    config = configparser.ConfigParser()
    return config.read('settings.config')

def get_address() -> [str, int]:
    config = get_configparser()
    return config['BASIC_PARAMS']['HOST'], int(config['BASIC_PARAMS']['PORT'])
