import json
import os
import _pickle as cPickle

def save(obj: object, prefixe: str = "output", format: str = "pickle"):

    assert format == 'json' or format == 'pickle', "File format is not supported"

    if format in ['json','pickle']:

        FILES = {'json': {'ext': 'json',
                          'mode': 'w',
                          'func': json.dump},
                'pickle': {'ext': 'pkl',
                          'mode': 'wb',
                          'func': cPickle.dump}}
        
        with open(file=f"{prefixe}.{FILES[format]['ext']}", mode=FILES[format]['mode']) as f:
            FILES[format]['func'](obj,f)

        assert os.path.isfile(f"{prefixe}.{FILES[format]['ext']}"), "File was not created"
    else:
        raise ValueError(f"Error: The file format {format} is not supported.")
    

def is_file(path: str) -> bool:
    return os.path.isfile(path)

def is_empty(path: str) -> bool:
    return os.path.getsize(path) == 0

def is_indexed(path: str) -> bool:
    return verify_file(path)

def verify_file(file: str):

    assert isinstance(file,str), "Values must be of string instance"

    if not is_file(file):

        raise FileNotFoundError(f"Error: The file '{file}' does not exist.")

    if is_empty(file):
            
        raise ValueError(f"Error: The file '{file}' is empty.")