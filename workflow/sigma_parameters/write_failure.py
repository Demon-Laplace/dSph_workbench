"""End the transfer protocol on a batch-level failure, preserving published data."""
import json,os,sys
from pathlib import Path
root=Path(__file__).resolve().parent
path=root/'transfer-final.json';temporary=root/'transfer-final.partial'
temporary.write_text(json.dumps(dict(terminal=True,final=True,expected_batches=1,completed_batches=1,
    batch_state='failed',exit_code=int(sys.argv[1]),published_data_preserved=True))+'\n')
os.replace(temporary,path)
