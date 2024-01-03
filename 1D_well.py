import numpy as np
import matplotlib.pyplot as plt

# Set the length of the box to 1 and number of electrons to 150
L = 1
num = 150

# Exact Energy Formula
def ex(n):
    return ((np.pi**2)/(24*L**2))*n*(n+2)*(n+1)

# The predicted energy has 2 contributions: the Von-Weizsacker (VW) energy, and the Harriman (H) correction term. Both must be numerically integrated.
# Harriman integrand:
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

# We integrate thes using the trapezoidal rule with resolution 1000.
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

# These formulae give the H and VW energy for a given electron number (n)

#Building arrays to plot:
electrons = np.arange(2,num+2,2)
exact= []
Har= []
V = []
pre = []
for i in range(len(electrons)):
    exact.append(ex(2*(i+1)))
    Har.append(H(2*(i+1)))
    V.append(VW(2*(i+1)))
    pre.append(VW(2*(i+1)) +H(2*(i+1)))

#Plotting
plt.scatter(electrons, Har, label = 'Harriman Correction')
plt.scatter(electrons, V, label = 'VW energy')
plt.scatter(electrons, exact, label = 'Exact energy')
plt.scatter(electrons, pre, label = 'Total predicted energy')
plt.ylabel('Energy')
plt.xlabel('Number of electrons')
plt.legend()
