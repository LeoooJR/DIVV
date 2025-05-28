import os
import pytest
import sys

tests_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(tests_dir)
if project_dir not in sys.path:
    sys.path.insert(0, project_dir)

import utils

@pytest.mark.utils
def test_runcmd():

    cmd: list[str] = ["echo", "is it working?"]
    
    out = utils.runcmd(cmd)

    assert not out.returncode

    assert out.stdout == b"is it working?\n"