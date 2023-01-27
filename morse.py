# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 20:38:09 2022

@author: Lakshya Singh
"""

#s-wave Schrodinger equation solution for V(r)$ = D(\exp(-2ar')-exp(-ar'))$
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

#inputting the value of all given constants
hbarc = 1973 #(eVA)
m = 940e6 #(eV/c^2)

#defining a function for potential
def pot(x):
    D = 0.755501 #eV
    a = 1.44 
    ro = 0.131349  # Angstrom
    r = (x-ro)/x
    return D*(np.exp(-2*r*a)-np.exp(-a*r))

#inputting range of x-axis and calculation of step size      
xo = 0.05
xn = 1.15
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

            
print("Kinetic energy matrix is given by: ",K) 
print()   
print("Potential energy matrix is given by: ",V) 

# calculating the hermition (H = K.E. + P.E.)
c = (-hbarc**2)/(2*m*h*h)

H = c*K + V
print()
print("Hamiltonian matrix obtained: ",H)
print()

#finding the eigenvalues and eigenvector for H(psi) = E(psi)
val,vect = eigh(H)

print("{:<5}{:<30}".format("State","  Eigenvalues"))
print("---------------------------")
for p in range(5):
      print("{:<5}{:<30}".format(p+1,val[p]))
          
print()
print("A few components of Eigenvectors are: ")
for p in range(5):
      print("For",p+1,"state: ",vect[0:5,p])

#plotting of eigenvectors which represent the state of particle (psi)
fig , ax = plt.subplots(figsize=(10,7))

plt.rcParams["figure.dpi"] = 1000

plt.subplot(2,2,1)
plt.title("x v/s $\psi(x)$")
plt.xlabel("x")
plt.ylabel("$\psi(x)$") 
plt.grid(True)
for i in range(0,3):
    plt.plot(x,(vect[:,i]),label = f"n = {i}")

#plt.plot(x[40:],pot(x[40:])/5, ls = "--", lw = 2, color = "k",label = "Potential")
plt.legend()

#plotting probabilities for states |(psi)|^2 
plt.subplot(2,2,2)
plt.title("x v/s $|\psi(x)|^2$")
plt.xlabel("x")
plt.ylabel("$|\psi(x)|^2$") 
plt.grid(True)
#plt.plot(x,pot(x))
for i in range(0,3):
    plt.plot(x,(vect[:,i])**2,label = f"n = {i}")

plt.subplot(2,2,3)
plt.plot(x[:],pot(x[:]), ls = "--", lw = 2, color = "k",label = "Potential")
plt.legend()


plt.suptitle("s-wave Schrodinger equation solution for V(r)$ = D(\exp(-2ar')-exp(-ar'))$",fontweight = "bold")
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

    

  





















