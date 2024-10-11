import logging
from config import LOGGING_LEVEL, LOGGING_FORMAT, LOG_FILE_PATH

def setup_logging():
    logging.basicConfig(
        level=getattr(logging, LOGGING_LEVEL),
        format=LOGGING_FORMAT,
        filename=LOG_FILE_PATH  # To log to a file
    )


def validate_data():
    pass

def i_o_operations():
    pass

