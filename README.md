# -*- coding: utf-8 -*-
"""
# Created on Tue Mar  3 14:57:41 2020

# @author: Adrian Kania, Krzysztof Sarapata, Michal Gucwa, and Anna Wójcik-Augustyn

# This software is released under the GNU General Public License
# Please cite:
# Adrian Kania, Krzysztof Sarapata, Michal Gucwa, and Anna Wójcik-Augustyn, "An algebraic approach to determine the coefficients of torsional potential - refinement of OPLS-AA force field for dimethyl phosphoric acid molecule"

Analytical solution of dihedral fitting

The script ‘K_fit.py’ calculates (fits) Fourier parameters for dihedral energy part (in force field equation).
As a input it takes the following data files. First, the energy of a molecule for different conformations is demanded.
It can be obtained using quantum chemistry methods. The second file should include energies for the same states 
but without dihedral part (parameters are set to zero). The last input files consist of dihedral angles of a molecule
which were considered (one dih file for every dihedral angle).

According to the formula in cited paper, the script calculates the estimated value of torsion part based on the given energy values.
The calculated parameters are selected in that way to minimize the square of the difference between corresponding energies
from the first and second file.
All torsion angles can be grouped according to their types.
For instance the torsion angle formed from OH.8-P.6-OS.10-CT atoms belongs to the same group as this one: CT.1-OS.5-P.6-OH.8.

The same set of four parameters K1, K2, K3, K4 applies to every dihedral of a given type.
The ‘K_fit.py’ script optionally expects two values ‘from’ and ‘to’ 
which set the range of the resulting parameters K. These optionally values 
and mandatory paths to data files the script gets from the standard input or text file with paths to QM energy,
MM energy and all dihedral angles files.
In the result, the script returns the set of four optimal Fourier or Ryckaert-Belleman (RB) parameters for each dihedral types 
and optionally set of four Fourier parameters from declared range for each dihedral types. 
The cited work includes the method of obtaining Ryckaert-Belleman parameters based on Fourier ones.


usage:

./K_fit.py from to < paths2files.txt  > coefficients.out

option:

from, to	    optionally two values ‘from’ and ‘to’ which set the range of the resulting parameters K.
				The calculated K parameters should be included inside this range. This option is only for Fourier parameters.

path2files.txt  includes paths to data files: 
                     path to QM energy file,
                     MM energy file,
                     *.dih files
                     and the amount of coefficients to calculate (3 or 4 value is allowed).

data files:

QM energy file	the energy of a molecule for different conformations e.g. quantum energy values from Gaussian 
				used as a base for making a fit. Values should have column layout.
MM energy file	the energy of a molecule for different conformations calculated with K's parameters set to 0. 
				Values should have column layout.


*.dih files 	with values of all dihedrals in molecule for different conformations. Values should have column layout. 
				The first line in each files containing the type of dihedral: order of the atoms forming that angle. 
				e.g. 'P.6-OS.10-CT.11-HC.14.dih' file in first line should have: 'P OS CT HC'.
				Ordering of values in all *.dih files is the same, so rows are related across all files.
coff            amount of coefficient to calculte. Default is 4, another possible one is 3

coffType        return the Ryckaert-Belleman coefficients instead of Fourier ones.

All QM, MM and *.dih files should include data in the same order which corresponds to the mutual relationship.
The path to those files are included in paths2files.txt

output:

In the result, the script returns the set of four optimal Fourier or RB parameters for each dihedral types 
and optionally set of four Fourier parameters from declared range for each dihedral types. 
The values are in JSON data interchange format.

"""
