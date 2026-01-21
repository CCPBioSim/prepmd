#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Change the modeller license key
"""

import importlib
from pathlib import Path
import sys

# note: don't do this
# i have to do it like this to make the conda build work

HELP_MSG =  "Usage: prep-license LICENSE-KEY (will replace current key)"

def entry_point():
    if len(sys.argv) != 2:
        print(HELP_MSG)
        sys.exit(0)
    if sys.argv[1] == "--help" or sys.argv[1] == "-h":
        print(HELP_MSG)
        sys.exit(0)
    key = sys.argv[1]
    modeller_init_path = Path(importlib.util.find_spec("modeller").origin)
    modeller_lib = modeller_init_path.parent.parent.absolute() / "modeller" / "config.py"
    with open(modeller_lib) as file:
        contents = file.readlines()
    splitline = contents[1].split("'")
    contents_1_new = splitline[0] + "'" + key + "'" + splitline[2]
    contents[1] = contents_1_new
    with open(modeller_lib, "w") as file:
        file.writelines(contents)
    print("Updated modeller license info.")

if __name__ == "__main__":

    entry_point()