"""
Add an offset to the obs_ids of a DL1 / DL2 file of ctapipe, in place: every table with an ``obs_id``
column (observation blocks, simulation runs, shower distributions, events). The simulations of several
zenith angles have the same CORSIKA runs (obs_ids), which ctapipe-merge cannot merge: mc_merge_train.bash
gives each zenith angle its own obs_ids (offset zenith * 10000000) in a copy of its training file.

    python offset_obs_ids.py <file> <offset>
"""
import sys

import numpy as np
import tables


def offset_obs_ids(path, offset):
    """Add ``offset`` to the column ``obs_id`` of all the tables of the HDF5 file ``path``; their number"""
    n_tables = 0
    with tables.open_file(path, mode="a") as h5file:
        for table in h5file.walk_nodes("/", classname="Table"):
            if "obs_id" not in table.colnames or table.nrows == 0:
                continue
            obs_ids = table.col("obs_id")
            if obs_ids.max() + offset > np.iinfo(obs_ids.dtype).max:
                raise OverflowError(f"{table._v_pathname}: obs_id + {offset} does not fit in {obs_ids.dtype}")
            table.modify_column(colname="obs_id", column=obs_ids + obs_ids.dtype.type(offset))
            n_tables += 1
    return n_tables


if __name__ == "__main__":
    path, offset = sys.argv[1], int(sys.argv[2])
    print(f"{path}: obs_id + {offset} in {offset_obs_ids(path, offset)} tables")
