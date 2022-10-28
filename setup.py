import setuptools
from setuptools import setup

from pybind11.setup_helpers import Pybind11Extension, build_ext


__version__ = "0.0.1"

cpp_args = ['-std=c++11', "-O3"]

ext_modules = [
    Pybind11Extension("_dfs_codes",
        ["cpp/graph_fast.cpp", "cpp/binding.cpp", "cpp/common.cpp"],# liz
        #["cpp/graph.cpp", "cpp/binding.cpp"],
        define_macros = [('VERSION_INFO', __version__)],
		extra_compile_args= cpp_args,
		extra_link_args= cpp_args)
]

setup(
    name="dfscode-wszola-wendler",
    version=__version__,
    author="Eliza Wszola & Chris Wendler",
    author_email="chris.wendler@inf.ethz.ch",
    description="Toolkit for the computation of (minimal) DFS codes.",
    long_description="",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    install_requires=['pybind11'],
    packages=setuptools.find_packages(where = 'src'),
    package_dir = {"":"src"}
)
