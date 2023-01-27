# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 17:20:25 2022

@author: Lakshya Singh
"""

#position and momentum commutator

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
    return sp.sin(x)**2

a = -i*hbar

def O(r):
    return sp.diff(funct(),r)

    
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
    if (B == O(r2) and A == O(r1)):
        ans1 = sp.diff(B,r1)
        ans2 = sp.diff(A,r2)
        return a*a*(ans1-ans2)
    
    elif (A == O(r1)):
        ans1 = sp.diff(B*funct(),r1)
        ans2 = B*A
        return a*(ans1-ans2)
    
    elif(B == O(r2)):
        ans2 = sp.diff(A*funct(),r2)
        ans1 = A*B
        return a*(ans1-ans2)
    else:
        return (A*B)*funct()-(B*A)*funct()    
    
A = O(x)
B = x
print("Value of Commutator = ",commute(A,B,x,x))  



























































'''
def Mcommute(r1,r2):
    #r1 diff wrt
    #A 2nd argument of commutator
    ans1 = sO.diff((r2*funct(x)),r1)
    ans2 = r2*(sp.diff(funct(x),r1))
    return a*(ans1-ans2)

def Mcommute(b,a):
    #b defines value of 2nd bracket in commutator bracket
    # a defines differtiation w.r.t   
    A =  Op1(b*f(x),a)
    B = b*Op1(f(x),x)
    return A-B
'''
   
 

"""
x,y,z = sp.symbols("x,y,z")

exp1 = sp.cos(x)
exp2 = sp.cos(-x)**2
print(exp1.subs(x,0))

f = sp.lambdify(x,exp1,np)
print(f(0))
print(sp.diff(exp1,x,8))
print(sp.integrate(exp1,(x,0,5)))
m = sp.integrate(exp2,(x,0,2))
print(m)
print(sp.limit(exp1,x,0))
p = exp1.series(x,0,4)
print(p)
print(sp.collect(exp2,x))
expr = (4*x**3 + 21*x**2 + 10*x + 12)/(x**4 + 5*x**3 + 5*x**2 + 4*x)
print(sp.apart(expr))

"""
      
      
      
      
      