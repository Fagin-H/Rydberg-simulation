# -*- coding: utf-8 -*-
"""
Created on Fri May 18 16:57:21 2018

@author: fhales
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt




#sigmap = sigmax() + sigmay()*1j
#sigmam = sigmax() - sigmay()*1j
I = qeye(2)
times = np.linspace(0.0, 81.0, 1000.0)
C3=1

psi1 = tensor(fock(2,1),fock(2,0),fock(2,0))
psi1 = psi1*psi1.dag()


H = 1/2*(C3/(1)*(tensor(sigmap(),sigmam(),I)+tensor(sigmam(),sigmap(),I))+C3/(1)*(tensor(I,sigmap(),sigmam())+tensor(I,sigmam(),sigmap()))+C3/(8)*(tensor(sigmap(),I,sigmam())+tensor(sigmam(),I,sigmap())))

psi0 = tensor(fock(2,1),fock(2,0),fock(2,0))

data = mcsolve(H,psi0,times,[],[psi1])

plt.plot(times, data.expect[0])