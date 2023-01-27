# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 14:43:24 2022

@author: Lakshya Singh
"""

import sympy as sp
#import numpy as np
from sympy.physics.quantum import Commutator, Operator
#from sympy.quantum.constants import h,pi,i

x,y,z,hbar,i,h = sp.symbols("x,y,z,hbar,i,h")


A = Operator("A")
B = Operator("B")
C = Operator("C")
D = Operator("D")


r1 = Commutator(A,B)
print(r1,"=",r1.doit())

r2 = Commutator(B,A)
print("[B,A] =",r2.doit(),"=",r2.expand(commutator = True))

r3 = Commutator(A,B+C+D)
print("[A,B+C+D] =",r3.doit(),"=",r3.expand(commutator = True))

r3 = Commutator(A,B-C-D)
print("[A,B-C-D] =",r3.doit(),"=",r3.expand(commutator = True))

print()
n = int(input("Enter the value of n for power property of Commutator: "))
r4 = Commutator(A,B**n)
r5 = Commutator(A**n,B)
print("[A,B^",n,"]","=",r4.doit(),"=",r4.expand(commutator = True))
print("[A^",n,"B] =",r5.doit(),"=",r5.expand(commutator = True))





