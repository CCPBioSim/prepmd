Python API
==========

``prepmd`` is also accessible via a Python API:

```python
from prepmd.add_modeller_license import setup_license
setup_license("<your license key>")  # only needs to be done once

from prepmd.prep import prep
prep("6xov", "6xov.cif", "working_dir")

from prepmd.run import run
run("6xov.cif", traj_out="traj.xtc")
```
