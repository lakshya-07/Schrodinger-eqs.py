# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 05:16:40 2022

@author: Lakshya Singh
"""

# angular momentum commutator

import sympy as sp
#import numpy as np
#from sympy.physics.quantum import Commutator, Dagger, Operator
#from sympy.quantum.constants import h,pi,i


x,y,z,hbar,i,h = sp.symbols("x,y,z,hbar,i,h")

def funct():
    '''
     defining psi.
    '''
    #return x**3 + 5*x
    #return sp.exp(x)
    return sp.sin(x)

a = -1j*hbar


def L(r):
    if (r == x):
        return y*a*sp.diff(funct(),z)-z*a*sp.diff(funct(),y)
    elif (r == y):
        return z*a*sp.diff(funct(),x)-x*a*sp.diff(funct(),z)
    elif (r == z):
        return x*a*sp.diff(funct(),y)-y*a*sp.diff(funct(),x)
    else:
        return "Your input is invalid. Choose r as x, y or z"
        
    
def commute(A,B,r1,r2):
    """
    Parameters
    ----------
    A : Term 1 of commutator bracket
        can be a variable or an operator
    B : Term 2 of commutator bracket
        can be a variable or an operator
    r1 : variable for which differentiation is done for operator A
        if A is variable then any value can be given to r1
    r2 : variable for which differentiation is done for operator B
        if B is variable then any value can be given to r2
    ----------
    """
    if (B == L(r2) and A == L(r1)):
        return hbar**2*(r1*sp.diff(funct(),r2)-r2*sp.diff(funct(),r1))
    
    elif(r1 == 2 or r2 == 2):
        return 0
    
    elif (A == L(r1)):
        if (r1 == x):
            ans1 = y*a*sp.diff(funct()*B,z)-z*a*sp.diff(B*funct(),y)
            ans2 = B*A
            return (ans1-ans2)
        elif (r1 == y):
            ans1 = z*a*sp.diff(B*funct(),x)-x*a*sp.diff(B*funct(),z)
            ans2 = B*A
            return (ans1-ans2)
        else:
            ans1 = x*a*sp.diff(B*funct(),y)-y*a*sp.diff(B*funct(),x)
            ans2 = B*A
            return (ans1-ans2)
    
    elif(B == L(r2)):
        if (r2 == x):
            ans1 = y*a*sp.diff(funct()*A,z)-z*a*sp.diff(A*funct(),y)
            ans2 = B*A
            return -(ans1-ans2)
        elif (r2 == y):
            ans1 = z*a*sp.diff(A*funct(),x)-x*a*sp.diff(A*funct(),z)
            ans2 = B*A
            return -(ans1-ans2)
        else:
            ans1 = x*a*sp.diff(B*funct(),y)-y*a*sp.diff(B*funct(),x)
            ans2 = B*A
            return -(ans1-ans2)
    else:
        return (A*B)*funct()-(B*A)*funct()    

    
A = z #sp.diff(f(x),y)
B = L(y)
print("Value of Commutator = ",commute(A,B,x,y))  





#print(-a*L(y))


#((L(x)**2*B - B*L(x)**2)+(L(y)**2*B - B*L(y)**2)+(L(z)**2*B - B*L(z)**2))
