import logging
from time import strftime, localtime
from enum import Enum


class LogLevel(Enum):
    """
    An enumeration representing different logging levels.
    Members:
        DEBUG (str): Represents a debugging log level, typically used for detailed diagnostic information.
        ERROR (str): Represents an error log level, used for logging error messages that indicate a problem in the program.
        INFO (str): Represents an informational log level, used for logging general system or application information.
        WARNING (str): Represents a warning log level, used for logging potentially harmful situations or messages indicating caution.
    """

    DEBUG = "debug"
    ERROR = "error"
    INFO = "info"
    WARNING = "warning"


def print_and_log(text: str, log_type: str):
    """
    Prints the given text with a timestamp and logs it with the specified log type.

    @param text: The text message to be printed and logged.
    @param log_type: The type of logging to be used. Valid values are 'debug', 'error', 'info', and 'warning'.

    This function first prints the provided text to the console, prefixed with a timestamp in the format '[YYYY-MM-DD HH:MM:SS]'. It then logs the same text using Python's logging module, categorized under the specified log type. The log type determines the severity level of the log message.

    """
    timestamp = strftime("[%Y-%m-%d %H:%M:%S]", localtime())
    print(f"{timestamp} {text}")

    if log_type == "debug":
        logging.debug(text)
    if log_type == "error":
        logging.error(text)
    if log_type == "info":
        logging.info(text)
    if log_type == "warning":
        logging.warning(text)
