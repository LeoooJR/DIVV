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

    checks = {files[0]: False,
              files[1]: False}

    for file in files:

        if is_file(file):

            if not is_empty(file):
                checks[file] = True
            else:
                pass
        else:
            pass

    return any(checks.values())