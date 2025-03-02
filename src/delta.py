import utils

def delta(params: object) -> int:

    print(params.vcfs)

    assert len(params.vcfs) == 2, "Two VCF files are required"

    if(utils.verify_files(files=params.vcfs)):
        pass
    else:
        return 1