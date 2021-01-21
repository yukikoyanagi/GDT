This file describes the files in GDT codebase.

Codes
-----
Each .py file contains documentation in the header.
Sample data is available in sample_data directory.

``compare_decoys_to_target.py`` computes P and S_n scores. ``const``, ``Calpha_max`` and ``Cbeta_max``  arguments are depreciated.

Make input file for ``least_squares.py`` using outputs from ``compare_decoys_to_target.py``. The input file should be a tab-delimited file, with each line contatining

``{dependent var} {independent var1} {independent var2} ...``

``gdtg_ts.py`` computes GDT-graph scores.

