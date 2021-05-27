from setuptools import setup

setup(
    name="GenTopo",
    version="0.1.0",
    description="Light-weight toolkit for topology manipulation",
    url="https://github.com/masrul/GenTopo",
    author="Masrul Huda",
    author_email="mmh568@msstate.edu",
    packages=["GenTopo"],
    py_modules=["Coord", "Graph", "GMXTopo", "PeriodicTable", "ImproperDihedral"],
    classifiers=[
        "Intended Audience :: Molecular Simulation",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.7",
    ],
)
