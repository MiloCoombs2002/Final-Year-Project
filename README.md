# Exact vs. Predicted Kinetic Energies in 1D Infinite Square Well
## Introduction
This repository showcases my final year project, focusing on the calculation of exact kinetic energies of electrons in a 1D infinite square well. The project explores the accuracy of predicted kinetic energies derived from the Harriman Orbitals proposed by John Harriman in 1981, employing orbital-free density functional theory.
## Key Features
- Calculation of exact kinetic energies in a 1D infinite square well.
- Utilization of Harriman Orbitals formula for predicting kinetic energies.
- Comparative analysis between predicted and exact kinetic energies to assess accuracy.

## Installation
1. Clone the repository to your local machine.
2. Ensure you have the required libraries installed (`numpy`, `matplotlib`).

## Usage
1. Run the provided code to calculate exact kinetic energies.
2. Apply the Harriman Orbitals formula to predict kinetic energies.
3. Compare and analyze the results using the provided tools.

## Code
#import libraries
import numpy as np
import matplotlib.pyplot as plt

#Set Length of box to 1
L = 1

#exact energy formula
def ex(n):
    return ((np.pi**2)/(24*L**2))*n*(n+2)*(n+1)

#The Predicted energy formual has two parts to to it; the Von-Weizsacker (VW) energy and the Harriman correction term that must both be numerically integrated using the trapezoidal rule.
#Harriman Integrand:
def Hi(n,x):
    a = 0
    for i in range(int(n/2)):
        a += np.sin(((i+1)*np.pi*x)/L)**2
    return a**3
#VW Integrand:
def VWi(n,x):
    a, b = 0, 0
    for i in range(int(n/2)):
        a += (i+1)*np.sin((i+1)*np.pi*x/L)*np.cos((i+1)*np.pi*x/L)
        b += np.sin((i+1)*np.pi*x/L)**2
    return (a**2)/(b+0.00000001)
#Integration with resolution of 1000. Th following forumale give H and VW energy for a given electron number:

def H(n):
    points = np.arange(0,1+1/res,1/res)
    s = Hi(n,points)
    s[0], s[-1] = s[0]/2, s[-1]/2
    integral = sum(s)/res
    return ((8*np.pi**2)/(3*L**3))*integral
def VW(n):
    points = np.arange(0,1+1/res,1/res)
    s = VWi(n,points)
    s[0], s[-1] = s[0]/2, s[-1]/2
    integral = sum(s)/res
    
    return (2*np.pi**2)*integral/(L**3)

#Building Arrays and Plotting Results
num = 150 #150 electrons

electrons = np.arange(2,num+2,2) # We only consider full orbitals (2 in each one)
exact= []
Har= []
V = []
pre = []
for i in range(len(electrons)):
    exact.append(ex(2*(i+1)))
    Har.append(H(2*(i+1)))
    V.append(VW(2*(i+1)))
    pre.append(VW(2*(i+1)) +H(2*(i+1)))

plt.scatter(electrons, Har, label = 'Harriman Correction')
plt.scatter(electrons, V, label = 'VW energy')
plt.scatter(electrons, exact, label = 'Exact energy')
plt.scatter(electrons, pre, label = 'Total predicted energy')
plt.xticks(ticks = electrons)
plt.ylabel('Energy')
plt.xlabel('Number of electrons')
plt.legend()

#Ratio Plot
ratio = np.array(pre)/np.array(exact)
plt.scatter(electrons,ratio,s=5)
plt.xlabel('Number of Electrons')
plt.ylabel('Ratio')
