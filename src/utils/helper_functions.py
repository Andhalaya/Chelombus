import logging
from config import LOGGING_LEVEL, LOGGING_FORMAT, LOG_FILE_PATH

def setup_logging():
    logging.basicConfig(
        level=getattr(logging, LOGGING_LEVEL),
        format=LOGGING_FORMAT,
        filename=LOG_FILE_PATH  # To log to a file
    )



def find_input_type(file_path):
    return file_path.split('.')[-1]

def validate_data():
    pass

def i_o_operations():
    pass

