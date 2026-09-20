import os
import platform
from pathlib import Path

from Cython.Build import cythonize
from setuptools import Extension, setup


def default_cplex_paths():
    if os.name == "nt":
        cplex_root = Path(
            r"C:\Program Files\IBM\ILOG\CPLEX_Studio2211\cplex"
        )
        include_dir = cplex_root / "include" / "ilcplex"
        library_dir = (
            cplex_root
            / "lib"
            / "x64_windows_msvc14"
            / "stat_mda"
        )
        libraries = ["cplex2211"]
        compile_args = []

    elif platform.system() == "Darwin":
        cplex_root = Path(
            "/Applications/CPLEX_Studio2211/cplex"
        )
        architecture = (
            "arm64_osx"
            if platform.machine() == "arm64"
            else "x86-64_osx"
        )
        include_dir = cplex_root / "include" / "ilcplex"
        library_dir = (
            cplex_root
            / "lib"
            / architecture
            / "static_pic"
        )
        libraries = ["cplex"]
        compile_args = ["-fPIC"]

    else:
        cplex_root = Path(
            "/opt/ibm/ILOG/CPLEX_Studio2211/cplex"
        )
        include_dir = cplex_root / "include" / "ilcplex"
        library_dir = (
            cplex_root
            / "lib"
            / "x86-64_linux"
            / "static_pic"
        )
        libraries = ["cplex"]
        compile_args = ["-fPIC"]

    return (
        include_dir,
        library_dir,
        libraries,
        compile_args,
    )


(
    default_include_dir,
    default_library_dir,
    libraries,
    compile_args,
) = default_cplex_paths()

cplex_include_dir = Path(
    os.environ.get(
        "CPLEX_INCLUDE_DIR",
        default_include_dir,
    )
)

cplex_library_dir = Path(
    os.environ.get(
        "CPLEX_LIBRARY",
        default_library_dir,
    )
)

if not cplex_include_dir.is_dir():
    raise FileNotFoundError(
        f"CPLEX include directory not found: "
        f"{cplex_include_dir}"
    )

if not cplex_library_dir.is_dir():
    raise FileNotFoundError(
        f"CPLEX library directory not found: "
        f"{cplex_library_dir}"
    )

print(f"CPLEX_INCLUDE_DIR={cplex_include_dir}")
print(f"CPLEX_LIBRARY={cplex_library_dir}")
print(f"CPLEX_LIBRARIES={libraries}")

extension = Extension(
    name="cplex_c_api_wrapper",
    sources=["cplex_c_api_wrapper.pyx"],
    include_dirs=[str(cplex_include_dir)],
    library_dirs=[str(cplex_library_dir)],
    libraries=libraries,
    extra_compile_args=compile_args,
)

setup(
    name="cplex_wrapper",
    ext_modules=cythonize(
        [extension],
        language_level=3,
    ),
)