# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 16:57:28 2022

@author: Lakshya Singh
"""
#s-wave Schrodinger equation solution for V(r) = (1/2)kr^2 + (1/3)br^3
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

#inputting the value of all given constants
hbarc = 197.3 #(MeV)
m = 940 #(MeV/c^2)

#defining a function for potential
bval = 11 # MeV fm^(-3) # given values are 0,10,30
def pot(x):
    k = 100 #Mev
    return (1/2)*k*x*x + (1/3)*bval*(x**3)

#inputting range of x-axis and calculation of step size      
xo = -5
xn = 5
n = 1001

x = np.linspace(xo,xn,n)
h = x[1]-x[0]

#forming Potential energy matrix
V = np.zeros((n,n))

for i in range(n):
    V[i][i] = pot(x[i])

#forming Kinetic energy matrix
K = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        if (i==j):
            K[i][i] = -2
        elif (abs(i-j)==1):
            K[i][j] = 1

            
print(f"Kinetic energy matrix for b = {bval} MeV fm^(-3) is given by: ",K) 
print()   
print(f"Potential energy matrix for b = {bval} MeV fm^(-3) is given by: ",V) 

# calculating the hermition (H = K.E. + P.E.)
c = (-hbarc**2)/(2*m*h*h)

H = c*K + V
print()
print(f"Hamiltonian matrix obtained for b = {bval} MeV fm^(-3): ",H)
print()

#finding the eigenvalues and eigenvector for H(psi) = E(psi)
val,vect = eigh(H)

print(f"Ground state energy for b = {bval}  MeV fm^(-3) is:",val[1],"MeV")
print(f"Ground state eigenvector for b = {bval} MeV fm^(-3) is: ",vect[:,1])


#plotting of eigenvectors which represent the state of particle (psi)
fig , ax = plt.subplots(figsize=(10,7))

plt.rcParams["figure.dpi"] = 1000

plt.subplot(2,2,1)
plt.title("x v/s $\psi(x)$")
plt.xlabel("x")
plt.ylabel("$\psi(x)$") 
plt.grid(True)

for i in range(1,4):
    plt.plot(x,(vect[:,i]),label = f"n = {i-1}")
    
#plt.plot(x,pot(x)*1e-4, ls = "--", lw = 2, color = "k",label = "Potential")
plt.legend()

#plotting probabilities for states |(psi)|^2 
plt.subplot(2,2,2)
plt.title("x v/s $|\psi(x)|^2$")
plt.xlabel("x")
plt.ylabel("$|\psi(x)|^2$") 
plt.grid(True)

for i in range(1,4):
    plt.plot(x,(vect[:,i])**2,label = f"n = {i-1}")
    
plt.legend()

plt.subplot(2,2,3)
plt.plot(x,pot(x), ls = "--", lw = 2, color = "k",label = "Potential")
plt.legend()


plt.suptitle(f"s-wave Schrodinger equation solution for V(r)$ = (1/2)kr^2 + (1/3)br^3$ for b = {bval} MeV fm^(-3)",fontweight = "bold")
plt.tight_layout(w_pad = 1, h_pad = -1)    
plt.show()













"""
    for j in range(0,len(x)):
        eigvect1.append(vect[j][i])
        eigvect2.append(vect[j][i]**2)
    
    #plt.subplot(1,2,1)
    plt.title("x v/s $\psi(x)$")
    plt.xlabel("x")
    plt.ylabel("$\psi(x)$") 
    plt.plot(x,pot(x))  
    plt.plot(x,eigvect1)

   
    plt.subplot(1,2,2)
    plt.title("x v/s $|\psi(x)|^2$")
    plt.xlabel("x")
    plt.ylabel("$|\psi(x)|^2$")    
    plt.plot(x,eigvect2)
    
  
    eigvect1 = []
    eigvect2 = []
    


plt.suptitle("Schrodinger wave solution for Harmonic Oscillator")
fig.legend(["n = 1","n = 2","n = 3"], loc = 'lower center')

plt.tight_layout(w_pad = 0.1, h_pad = -1)


"""

    

  





















