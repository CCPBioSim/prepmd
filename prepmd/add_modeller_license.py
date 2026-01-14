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

def entry_point():
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
    if len(sys.argv) != 2:
        print("Usage: prep-license LICENSE-KEY")
        print("Note: this will replace whatever your current license key is.")
        sys.exit(0)
    entry_point()