from setuptools import setup, Extension
from Cython.Build import cythonize
import os

CPLEX_INCLUDE_DIR = os.environ.get("CPLEX_INCLUDE_DIR", "/opt/ibm/ILOG/CPLEX_Studio/cplex/include")
CPLEX_LIBRARY = os.environ.get("CPLEX_LIBRARY", "/opt/ibm/ILOG/CPLEX_Studio/cplex/lib/x86-64_linux/static_pic")

print("CPLEX_LIBRARY=", CPLEX_LIBRARY)
print("CPLEX_INCLUDE_DIR=", CPLEX_LIBRARY)

ext_modules = [
    Extension(
        name="cplex_c_api_wrapper",
        sources=["cplex_c_api_wrapper.pyx"],
        include_dirs=[CPLEX_INCLUDE_DIR],
        library_dirs=[CPLEX_LIBRARY],
        libraries=["cplex"],
        extra_compile_args=["-fPIC"],
    )
]

setup(
    name="cplex_wrapper",
    ext_modules=cythonize(ext_modules, language_level=3),
)
