# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 21:15:13 2022

@author: Lakshya Singh
"""

#Schrodinger wave solutions for Box Potential
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

#inputting the value of all given constants
hbar = 1 #1.05457182e-34
m = 1 #9.1093837e-31 

#inputting range of x-axis and calculation of step size
xo = 0 #-0.05
xn = 10 #0.05

n = 1001

x = np.linspace(xo,xn,n)

h = (x[1]-x[0]) 

#forming Potential energy matrix (which is zero for this case)
V = np.zeros((n,n))

for i in range(n):
    for j in range(n):
        if (i==j):
            w = 0
            V[i][i] = 0.5*m*w*w*(x[i]**2)
      
#forming Kinetic energy matrix
K = np.zeros((n,n))

for i in range(n):
    for j in range(n):
        if (i==j):
            K[i][i] = -2
        elif (abs(i-j)==1):
            K[i][j] = 1

            
print("Kinetic energy matrix is given by: ",K) 
print()   
print("Potential energy matrix is given by: ",V) 


# calculating the hermition (H = K.E. + P.E.)
c = (-hbar**2)/(2*m*h*h)
H = c*K + V
print()
print("Hamiltonian matrix obtained: ",H)
print()

#finding the eigenvalues and eigenvector for H(psi) = E(psi)
val,vect = eigh(H)

print("{:<5}{:<30}".format("State","  Eigenvalues"))
print("---------------------------")
for p in range(10):
      print("{:<5}{:<30}".format(p+1,val[p]))
          
print()
print("A few elements of Eigenvectors are: ")
for p in range(10):
      print("For",p,"state: ",vect[0:5,p])
      

#plotting of eigenvectors which represent the state of particle (psi)
fig , ax = plt.subplots(figsize=(7,7))
plt.rcParams["figure.dpi"] = 1000


plt.subplot(1,2,1)
plt.title("x v/s $\psi(x)$")
plt.xlabel("x")
plt.ylabel("$\psi(x)$") 
for i in range(0,6):
    plt.plot(x,(vect[:,i])+ i/5)

for i in range(0,11,2):
    plt.axhline(y = i/10,linestyle ="--",color = "grey")
    plt.text(xn-h ,i/10, f"n = {int(i/2)}",fontweight = "bold")


#plotting probabilities for states |(psi)|^2 
plt.subplot(1,2,2)
plt.title("x v/s $|\psi(x)|^2$")
plt.xlabel("x")
plt.ylabel("$|\psi(x)|^2$") 

for i in range(0,6):
    plt.plot(x,(vect[:,i])**2 + i/200)
    plt.text(xn-h ,i/200, f"n = {i}",fontweight = "bold")

   
for i in range(0,30,5):
    plt.axhline(y = i/1000,linestyle ="--",label = f"n = {i}",color = "grey")


plt.suptitle("Schrodinger wave solutions for Box Potential",fontweight = "bold")
plt.tight_layout(w_pad = 1, h_pad = -1)    
plt.show()








