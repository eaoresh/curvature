from setuptools import setup

from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "my_module",            # python module name goes here
        ["src/main.cpp"],       # relative path to cpp file
        cxx_std=17,             # required cxx_std version
    ),
]

setup(
    name="my_module",           # python module name goes here too
    ext_modules=ext_modules,
)
