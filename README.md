# lmp2cg

lmp2cg is a python tool for converting LAMMPS data files with explicit bonds to a coarse grained representation. For this purpose a template file is used. lmp2cg then searches in the data file for corresponding stuctures and outputs an coarse grained data file. Additionally dump files can converted. At the moment periodic boundary conditions are assumed.

Limitations:
- Unfortunately at the moment only unwrapped coordinates in both the data as well as the dump file are necessary to get correct coordinates.
- at the moment only atom_style "full" is supported.

This code is far from being finished, well documented or ready for release. It could be called beta (if at all).

Typical usage:
python3 lmp2cg --data lmp.data --template template.dat [--trajectory lmp.dump]

Python modules used:
sys, os, datetime, numpy, NetworkX
