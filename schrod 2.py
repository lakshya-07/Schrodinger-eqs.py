# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 14:25:32 2022

@author: Lakshya Singh
"""
#Schrodinger wave solutions for Harmonic Oscillator
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

#inputting the value of all given constants
hbar = 1.05457182e-34
m = 9.1093837e-31

#defining a function for potential
def pot(x):
    w = 1
    return 0.5*m*w*w*(x**2)

#inputting range of x-axis and calculation of step size  
xo = -0.05
xn = 0.05

n = 1001

x = np.linspace(xo,xn,n)

h = x[1]-x[0]

#forming Potential energy matrix
V = np.zeros((n,n))

for i in range(n):
    w = 1
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
print("A few components of Eigenvectors are: ")
for p in range(10):
      print("For",p,"state: ",vect[0:5,p])

#plotting of eigenvectors which represent the state of particle (psi)
fig , ax = plt.subplots(figsize=(7,7))
plt.rcParams["figure.dpi"] = 1000

plt.subplot(1,2,1)
plt.title("x v/s $\psi(x)$")
plt.xlabel("x")
plt.ylabel("$\psi(x)$") 
plt.plot(x,pot(x)*1e33)
for i in range(0,6):
    plt.plot(x,(vect[:,i])+ i/5)
    plt.text(xn-h ,i/5, f"n = {i}",fontweight = "bold")

for i in range(0,11,2):
    plt.axhline(y = i/10,linestyle ="--",color = "grey")    

#plotting probabilities for states |(psi)|^2 
plt.subplot(1,2,2)
plt.title("x v/s $|\psi(x)|^2$")
plt.xlabel("x")
plt.ylabel("$|\psi(x)|^2$") 
plt.plot(x,pot(x)*1e32)
for i in range(0,6):
    plt.plot(x,(vect[:,i])**2 + i/50)
    plt.text(xn-h ,i/50, f"n = {i}",fontweight = "bold")


for i in range(0,11,2):
    plt.axhline(y = i/100,linestyle ="--",label = f"n = {i}",color = "grey")

plt.suptitle("Schrodinger wave solutions for Harmonic Oscillator",fontweight = "bold")
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





K = np.diag([-2]*n,0)

for i in range(1,2):
    K += np.diag([1]*(n-i),i)
    K += np.diag([1]*(n-i),-i)

"""

    

  





















