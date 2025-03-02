import os

def is_file(path: str) -> bool:
    return os.path.isfile(path)

def is_empty(path: str) -> bool:
    return os.path.getsize(path) == 0

def is_indexed(path: str) -> bool:
    pass

def verify_files(files: list[str]):

    assert len(files) == 2, "Two files are required"
    assert isinstance(files[0],str) and isinstance(files[1],str), "Values must be of string instance"

    for file in files:

        if not is_file(file):

            raise FileNotFoundError(f"Error: The file '{file}' does not exist.")

        if is_empty(file):
            
            raise ValueError(f"Error: The file '{file}' is empty.")