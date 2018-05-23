# -*- coding: utf-8 -*-
"""
Created on Fri May 18 16:57:21 2018
@author: fhales
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randn

I = qeye(3)
times = np.linspace(0.0, 21, 1000.0)
C3 = 7965
kB = 1.38e-23
m = 1.42e-25
omega = 0.09
Temp = 5e-5

gamma = 1
gammaup = 1/101
gammadown = 1/135

sx = Qobj([[0,1,0],
           [1,0,0],
           [0,0,0]])

sy = Qobj([[0,-1j,0],
           [1j,0,0],
           [0,0,0]])

sigmap = sx + sy*1j
sigmam = sx - sy*1j

qubits = np.array([[0,0,0],[20,0,0],[40,0,0]])

opts = Odeoptions(nsteps=500000)

def makematrix(qubits_co):
    intmatrix = [[np.linalg.norm(qubits_co[i]-qubits_co[j])**3 for j in range(len(qubits_co))] for i in range(len(qubits_co))]
    return intmatrix

def makehamcoeff(qubits,q1,q2,T):
    r_rms = (T*kB/(m*omega**2))**0.5
    v_rms = (T*kB/(m))**0.5
    q1pos = qubits[q1] + 3**-0.5*r_rms*randn(3)
    q2pos = qubits[q2] + 3**-0.5*r_rms*randn(3)
    def H_coeff(t,args):
        return 0.5*C3*np.linalg.norm(q1pos + 3**-0.5*v_rms*randn(3)*t - q2pos - 3**-0.5*v_rms*randn(3)*t)**-3
    return H_coeff

def makeinputoutput(atom_number):
    tempstatein = fock(3,1)
    for i in range(atom_number-1):
        tempstatein = tensor(tempstatein,fock(3,0))
    tempstateout = tempstatein*tempstatein.dag()
    return tempstatein, tempstateout

def makesig(i,j,atoms):
    temp1 = [I]*atoms
    temp2 = [I]*atoms
    temp1[i] = sigmap
    temp1[j] = sigmam
    temp2[i] = sigmam
    temp2[j] = sigmap
    temp = tensor(temp1)+tensor(temp2)
    return temp

def makecollapse(i,qubits):
    c_ops = []
    temp = [I]*len(qubits)
    temp[i] = (gamma+gammaup)**0.5*fock(3,2)*fock(3,1).dag()
    tens1 = tensor(temp)
    c_ops.append(tens1)
    temp[i] = gammadown**0.5*fock(3,2)*fock(3,0).dag()
    tens2 = tensor(temp)
    c_ops.append(tens2)
    return c_ops

def makeham(intmatrix):
    atoms = len(intmatrix)
    components =  []
    for i in range(atoms):
        for j in range(atoms):
            if j < i:
                components.append(1/intmatrix[i][j]*makesig(i,j,atoms))
    ham = 0.5*C3*sum(components)
    return ham

def makehamt(qubits,T):
    return [tensor([0*I]*len(qubits))] + [[makesig(i,j,len(qubits)),makehamcoeff(qubits,i,j,T)] 
                   for i in range(len(qubits)) for j in range(len(qubits)) if j<i]
    

def makeplot(H,timesteps,instate,basis,qubits,collapse=None):
    c_ops = []
    if collapse == True:
        for i in range(len(qubits)):
            c_ops += makecollapse(i,qubits)
    else:
        pass
    data = mesolve(H,instate,timesteps,c_ops,basis,options=opts)
    plt.plot(timesteps, data.expect[0])

def doall(qubits,times,collapse=None):
    int_matrix = makematrix(qubits)
    input_state, output_basis = makeinputoutput(len(qubits))
    H = makeham(int_matrix)
    makeplot(H,times,input_state,output_basis,qubits,collapse)
    
def doallt(qubits,times,T,collapse=None):
    input_state, output_basis = makeinputoutput(len(qubits))
    H = makehamt(qubits,T)
    makeplot(H,times,input_state,output_basis,qubits,collapse)