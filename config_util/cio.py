import configparser
import os

config = configparser.ConfigParser()
path_current_directory = os.path.dirname(__file__)
path_config_file = os.path.join(path_current_directory, 'vc.config')
config.read(path_config_file)


# Basic params
def get_address() -> (str, int):
    """
    Gets field from vc.config_util
    @return: tuple consisting of host-IP as string and port as int
    """
    return config['BASIC_PARAMS']['HOST'], int(config['BASIC_PARAMS']['PORT'])


def get_queue_size() -> int:
    """
    Gets field from vc.config_util
    @return: queue size as int
    """
    return int(config['BASIC_PARAMS']['QUEUE_SIZE'])


def get_min_queue_size() -> int:
    """
    Gets field from vc.config_util
    @return: min queue size as int
    """
    return int(config['BASIC_PARAMS']['MIN_QUEUE_SIZE'])


def get_max_queue_size() -> int:
    """
    Gets field from vc.config_util
    @return: max queue size as int
    """
    return int(config['BASIC_PARAMS']['MAX_QUEUE_SIZE'])


def get_temp_dir() -> str:
    """
    Gets field from vc.config_util
    @return: path for temp files
    """
    return config['BASIC_PARAMS']['TEMP_DIR']


def get_temp_file_extension() -> str:
    """
    Gets field from vc.config_util
    @return: extension for temp files
    """
    return config['BASIC_PARAMS']['TEMP_FILE_EXTENSION']


# Variant Caller Params

def get_reference() -> str:
    """
    Gets field from vc.config_util
    @return: path to FASTA-reference file
    """
    return config['VARIANT_CALLER_PARAMS']['reference']


def get_min_evidence_depth() -> int:
    """
    Gets field from vc.config_util
    @return: minimal evidence depth as int
    """
    return int(config['VARIANT_CALLER_PARAMS']['minEvidenceDepth'])


def get_min_evidence_ratio() -> float:
    """
    Gets field from vc.config_util
    @return: minimal evidence ration as int
    """
    return float(config['VARIANT_CALLER_PARAMS']['minEvidenceRatio'])


def get_max_variants() -> int:
    """
    Gets field from vc.config_util
    @return: max. count of variants as int
    """
    return int(config['VARIANT_CALLER_PARAMS']['maxVariants'])


def get_min_total_depth() -> int:
    """
    Gets field from vc.config_util
    @return: minimal total depth as int
    """
    return int(config['VARIANT_CALLER_PARAMS']['minTotalDepth'])


def get_min_mapping_quality() -> int:
    """
    Gets field from vc.config_util
    @return:
    """
    return int(config['VARIANT_CALLER_PARAMS']['minMappingQuality'])


def get_min_base_quality() -> int:
    """
    Gets field from vc.config_util
    @return: minimal base quality as int
    """
    return int(config['VARIANT_CALLER_PARAMS']['minBaseQuality'])


# Watcher Params

def get_watcher_interval() -> int:
    """
    Gets field from vc.config_util
    @return: intervals in which dir will be checked in seconds
    """
    return int(config['WATCHER_PARAMS']['WATCHER_INTERVAL'])


def get_watch_recursively() -> bool:
    """
    Gets field from vc.config_util
    @return: boolean for recursive watching
    """
    return bool(config['WATCHER_PARAMS']['WATCH_RECURSIVELY'])


def get_supported_extensions() -> list[str]:
    """
    Gets field from vc.config_util
    @return: list of supported extensions by watcher script
    """
    return config['WATCHER_PARAMS']['SUPPORTED_EXTENSIONS'].split(',')
