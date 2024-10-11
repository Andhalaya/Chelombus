import logging
from config import LOGGING_LEVEL, LOGGING_FORMAT, LOG_FILE_PATH

def setup_logging():
    logging.basicConfig(
        level=getattr(logging, LOGGING_LEVEL),
        format=LOGGING_FORMAT,
        filename=LOG_FILE_PATH  # To log to a file
    )



def find_input_type(self, file_path):
    if file_path.endswith('csv'):
        return True
    # TODO: Add support for txt  
    # elif file_path.split('.')[-1] == 'txt':
    #     return False 
    else: 
        raise ValueError('Unsupported input file. Only .csv files are supported')
    

def validate_data():
    pass

def i_o_operations():
    pass

