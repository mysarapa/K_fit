#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
# Created on Tue Mar  3 14:57:41 2020

# @author: Adrian Kania, Anna Wójcik-Augustyn, Krzysztof Sarapata, …

# This software is released under the GNU General Public License
# Please cite:
# Kania A., Wójcik-Augustyn A., Sarapata K., Gucwa M., 
# An algebraic approach to determine the torsional coefficients - refinement of OPLS-AA force field for dimethyl phosphate molecule,
# 2020

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
In the result, the script returns the set of four optimal Fourier or Ryckaert-Belleman parameters for each dihedral types 
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

In the result, the script returns the set of four optimal Fourier parameters for each dihedral types 
and optionally set of four Fourier parameters from declared range for each dihedral types. 
The values are in JSON data interchange format.

"""


import os
import sys
import string
import pandas as pd
import numpy as np
import itertools
from numpy.linalg import inv
import json
from scipy.optimize import lsq_linear


print ("Analytic method of dihedral fitting")
print ("Adrian Kania")
print ("Departemnt of Computational Biophysics and Bioinformatics")
print ("Jagiellonian University")
print ("2020")
print ("This software is released under the GNU General Public License")
print ("Please cite:")
print ("   ")

def forDerivative(M,x): 
    '''for derivative calculation'''
    for i in range(len(x)):
        M[i]+=[x[i]*t for t in x]
    return M

def fourier2RM (k1,k2,k3,k4=0):
    '''Return the Ryckaert-Belleman coefficients for the given Fourier ones'''
    return k1+k3+2*k4, 3*k3-k1, 2*k2-8*k4, -4*k3, 8*k4, 0

# setting the range of calculated parameters: from to (If we want to have parmeters Kij in a certain range)

if len(sys.argv) == 3:
	fr = int(sys.argv[1])
	to = int(sys.argv[2])

print("What is the path of the qm energy file?")
qmFileName = sys.stdin.readline().strip()
qmFile = pd.read_csv(qmFileName, sep="\s+", header=None)

#Energy QM (Gaussian)
qm = qmFile[qmFile.columns[0]]
print("  qm energy file = {} was read".format(qmFileName))

#This file should come from fit_wrapper (gromax) and it contains four columns: qm energy, free energy (gromax) for zero dihedral part, energy after fitting 
print('What is the path of the mm energy file (second column)?')
mmFileName = sys.stdin.readline().strip()

# Energy from Gromax where torsion part is set to zero (fit_wrapper)

mmFile = pd.read_csv(mmFileName, sep="\s+", header=None)
zero_rb = mmFile[mmFile.columns[1]]
print("  mm energy file = {} was read".format(mmFileName))

print('What is the path to *.dih files?')
dihDir = sys.stdin.readline().strip()


# read *.dih files name from dir
files = [f for f in os.listdir(dihDir) if f.endswith('.dih')]
print('{} files *.dih were founded'.format(len(files)))


# type of coefficient
print('What type of coefficients Fourier/Ryckaert - Bellemans return?')
coffType = sys.stdin.readline().strip()
coffType = 'RB' if coffType == 'RB' else 'Fourier'
print('{} coefficients will be calculated'.format(coffType))


# read amount of output coefficients
print('What is the amount of coefficients?')
coff = sys.stdin.readline().strip()
coff = 3 if coff == '3' else 4
print('{} coefficients will be calculated'.format(coff))


diF = {k:open(dihDir + '/' + k).readline().strip() for k in files}

# recognition of angle types
ir = {}
for f in diF:#dih.columns)):
    if diF[f] not in ir and ' '.join(diF[f].split()[::-1]) not in ir:
        ir[diF[f]]=[f]
    else:
        if diF[f] in ir:
            ir[diF[f]].append(f)
        else:
            ir[' '.join(diF[f].split()[::-1])].append(f)
how_many=[]
files_orders = []
type_orders = []

# ordering files, types of angles
for k,v in ir.items():
    how_many.append(len(v))
    files_orders += v
    type_orders.append(k)


# data frame of all angles basis on files *.dih

fs = [pd.read_csv(dihDir + '/' + f) for f in files_orders]

dih = pd.concat(fs, axis=1)


n = 1
components = 0

B_t=(qm-zero_rb) #torsion component


sum_coefficients_K=[[] for x in range(len(how_many))]

WSP=[]
B=np.zeros((coff*len(how_many)))
b1=0
M=np.zeros((coff*len(how_many), coff*len(how_many)))
ss=0
i=0
cc=[]
#Kt=[]
#sup =[]
for index, row in dih.iterrows():
    for kat in row:
        for cos in range(coff):
            i+=1
            cc.append(kat)
            WSP.append((1+np.cos(n*(kat) * np.pi / 180 )))
            n+=1
            components += 1
            if n > coff:
                n = 1
            if components > len(dih.columns) * coff - 1:
                components = 0
                Kt =[WSP[x::coff] for x in range(coff)]
                tt = 0

                for j in range(len(how_many)):
                    for x in range(coff):
                        sum_coefficients_K[j].append(sum(Kt[x][tt:(tt + how_many[j])]))

                    tt = tt + how_many[j]


                sum_coefficients_K = list(itertools.chain(*sum_coefficients_K))  # flat na liscie
    
                M=forDerivative(M,sum_coefficients_K)

                for j in range(len(sum_coefficients_K)):
                    B[j]+=B_t[ss]*sum_coefficients_K[j]

                ss+=1

                WSP = []
                sum_coefficients_K=[[] for x in range(len(how_many))]

M=np.array(M)
B=np.array(B)
K=np.dot(inv(M), B)

result = {type_orders[i]:list(K)[coff*i:coff*i+coff] for i in range(len(type_orders))}

cT = {'RB': 'Ryckaert - Bellemans', 'Fourier': 'Fourier'}


print('\nThe best result: {} {} parametres for each type of dihedral angle:'.format(coff, cT[coffType]))

if coffType == 'Fourier':
    print(json.dumps(result, indent=1))
else:
    print(json.dumps({key: list(fourier2RM(*result[key])) for key in result}, indent=1))

if 'fr' in locals():
    x1 = [fr] * coff * len(how_many)
    x2 = [to] * coff * len(how_many)
    ww=lsq_linear(M,B,bounds=(x1, x2)).x
    print('\n')
    resultWW = {type_orders[i]:list(ww)[coff*i:coff*i+coff] for i in range(len(type_orders))}
    print('\nThe best result {} {} parametres in range ({}:{}) for each type of dihedral angle:'.format(coff, cT[coffType], fr, to))
    if coffType == 'Fourier':
        print(json.dumps(resultWW, indent=1))
    else:
        print(json.dumps({key: list(fourier2RM(*resultWW[key])) for key in result}, indent=1))

