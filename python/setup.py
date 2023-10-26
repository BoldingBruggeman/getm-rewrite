import os
import sys
from setuptools import setup, find_packages

FABM_BASE = os.path.join(os.path.abspath(os.path.dirname(__file__)), "../extern/fabm")
sys.path.insert(0, os.path.join(FABM_BASE, "src/drivers/python"))

from build_cmake import CMakeExtension, bdist_wheel, CMakeBuild

setup(
    name="pygetm",
    version="0.1.1",
    author="Bolding-Bruggeman ApS",
    author_email="jorn@bolding-bruggeman.com",
    license="GPL",
    packages=find_packages(include=["pygetm*"]),
    package_data={
        "pygetm": ["*.so", "*.dll", "*.dylib", "*.pyd"],
        # "pygetm.pyfabm": ["*.so", "*.dll", "*.dylib", "*.pyd"],
    },
    ext_modules=[
        CMakeExtension(
            "pygetm.fabm",
            "-DPYFABM_DEFINITIONS=_FABM_DIMENSION_COUNT_=3;_FABM_DEPTH_DIMENSION_INDEX_=3;_FABM_MASK_TYPE_=integer;_FABM_MASKED_VALUE_=0;_FABM_HORIZONTAL_MASK_;_FABM_VERTICAL_BOTTOM_TO_SURFACE_",
        ),
    ],
    zip_safe=False,
    entry_points={
        "console_scripts": [
            "pygetm-subdiv = pygetm.subdiv:main",
            "pygetm-test-scaling = pygetm.parallel:test_scaling_command",
            "pygetm-compare-nc = pygetm.util.compare_nc:compare_command",
        ],
    },
    cmdclass={"bdist_wheel": bdist_wheel, "build_ext": CMakeBuild},
)
