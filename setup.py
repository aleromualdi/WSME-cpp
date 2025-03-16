import pybind11
from setuptools import Extension, setup

ext_modules = [
    Extension(
        "wsme-cpp",
        ["cpp/wrapper.cpp", "cpp/partition_function.cpp"],
        include_dirs=[
            pybind11.get_include(),
            pybind11.get_include(True),
        ],
        language="c++",
        extra_compile_args=["-std=c++11"],
    ),
]

setup(
    name="wsme",
    ext_modules=ext_modules,
    zip_safe=False,
)
