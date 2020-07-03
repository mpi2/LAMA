"""
Get a checksum for the contents of a config file. This can be used to make sure unforseen errors creeep into config
"""
import hashlib
import json
from typing import Dict

def md5(data: Dict) -> str:
    return hashlib.md5(json.dumps(data, sort_keys=True)).hexdigest()
