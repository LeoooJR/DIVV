import utils

def delta(params: object) -> int:

    print(params.vcfs)

    assert len(params.vcfs) == 2, "Two VCF files are required"

    utils.verify_files(files=params.vcfs)

    if(params.indexes):

        assert len(params.indexes) == 2; "Two indexes files are required"
        
        utils.verify_files(files=params.indexes)
            
    else:
        pass
    
    return 1