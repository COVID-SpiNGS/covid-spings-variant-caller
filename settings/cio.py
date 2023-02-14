import configparser
import os

config = configparser.ConfigParser()
path_current_directory = os.path.dirname(__file__)
path_config_file = os.path.join(path_current_directory, 'vc_settings.config')
config.read(path_config_file)


# Basic params
def get_address() -> (str, int):
    return config['BASIC_PARAMS']['HOST'], int(config['BASIC_PARAMS']['PORT'])


def get_queue_size() -> int:
    return int(config['BASIC_PARAMS']['QUEUE_SIZE'])


def get_temp_dir() -> str:
    return config['BASIC_PARAMS']['TEMP_DIR']


def get_temp_file_extension() -> str:
    return config['BASIC_PARAMS']['TEMP_FILE_EXTENSION']


# Variant Caller Params

def get_reference() -> str:
    return config['VARIANT_CALLER_PARAMS']['reference']


def get_min_evidence_depth() -> int:
    return int(config['VARIANT_CALLER_PARAMS']['minEvidenceDepth'])


def get_min_evidence_ratio() -> float:
    return float(config['VARIANT_CALLER_PARAMS']['minEvidenceRatio'])


def get_max_variants() -> int:
    return int(config['VARIANT_CALLER_PARAMS']['maxVariants'])


def get_min_total_depth() -> int:
    return int(config['VARIANT_CALLER_PARAMS']['minTotalDepth'])


def get_min_mapping_quality() -> int:
    return int(config['VARIANT_CALLER_PARAMS']['minMappingQuality'])


def get_min_base_quality() -> int:
    return int(config['VARIANT_CALLER_PARAMS']['minBaseQuality'])


# Watcher Params

def get_watcher_interval() -> int:
    return int(config['WATCHER_PARAMS']['WATCHER_INTERVAL'])


def get_watch_recursively() -> bool:
    return bool(config['WATCHER_PARAMS']['WATCH_RECURSIVELY'])


def get_supported_extensions() -> list[str]:
    return config['WATCHER_PARAMS']['SUPPORTED_EXTENSIONS'].split(',')
