import os
from pathlib import Path


_CPLEX_DLL_DIRECTORY = None

if os.name == "nt":
    default_cplex_bin = Path(
        r"C:\Program Files\IBM\ILOG\CPLEX_Studio2211"
        r"\cplex\bin\x64_win64"
    )

    cplex_bin = Path(
        os.environ.get(
            "CPLEX_BIN_DIR",
            str(default_cplex_bin),
        )
    )

    if not cplex_bin.is_dir():
        raise FileNotFoundError(
            f"CPLEX DLL directory not found: {cplex_bin}"
        )

    # Keep the handle alive for the lifetime of this package.
    _CPLEX_DLL_DIRECTORY = os.add_dll_directory(
        str(cplex_bin)
    )

from .cplex_c_api_wrapper import *