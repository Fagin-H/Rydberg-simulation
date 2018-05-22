# -*- coding: utf-8 -*-
"""
Created on Fri May 18 16:57:21 2018
@author: fhales
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

I = qeye(2)
times = np.linspace(0.0, 1.0, 100.0)
C3 = 1


initqubits = [0,0,0,1,1,0,1,1,0,2]

#outputs = [tensor(fock(2,0),fock(2,0))*tensor(fock(2,0),fock(2,0)).dag(),tensor(fock(2,0),fock(2,1))*tensor(fock(2,0),fock(2,1)).dag(),tensor(fock(2,1),fock(2,1))*tensor(fock(2,1),fock(2,1)).dag(),tensor(fock(2,1),fock(2,0))*tensor(fock(2,1),fock(2,0)).dag()]
outputs = [fock(2,1),fock(2,0)]

def makematrix(qubits_co):
    intmatrix = []
    for i in range(len(qubits_co)):
        temprow = []
        for j in range(len(qubits_co)):
            temprow.append(((qubits_co[i][0]-qubits_co[j][0])**2+(qubits_co[i][1]-qubits_co[j][1])**2)**1.5)
        intmatrix.append(temprow)
    return intmatrix
        

def makesig(i,j,atoms):
    temp1 = [I]*atoms
    temp2 = [I]*atoms
    temp1[i] = sigmap()
    temp1[j] = sigmam()
    temp2[i] = sigmam()
    temp2[j] = sigmap()
    temp = tensor(temp1)+tensor(temp2)
    return temp
    

def makeham(intmatrix):
    atoms = len(intmatrix)
    components =  []
    for i in range(atoms):
        for j in range(atoms):
            if j < i:
                components.append(1/intmatrix[i][j]*makesig(i,j,atoms))
    ham = C3*sum(components)
    return ham


def makeplot(H,timesteps,instate,basis):
    data = mcsolve(H,instate,timesteps,[],basis)
    plt.plot(timesteps, data.expect[0])

def getdata(qubits,times,input_state):
    input_state = tensor(input_state,fock(2,0),fock(2,0),fock(2,1),fock(2,1))
    int_matrix = makematrix(qubits)
    H = makeham(int_matrix)
    data = mcsolve(H,input_state,times,[],[]).states[-1].ptrace([0])
    return data

def makeallstates(n):
    states = []
    for i in range(n):
        binn = bin(i)[2:]
        while len(binn) != len(bin(n)[2:])-1:
            binn = '0' + binn
        tempstate = [fock(2,int(binn[j])) for j in range(len(binn))]
        states.append(tensor(tempstate))
    return states

def calfil(qubits):
    #qubits = [[qubits[n],0] for n in range(len(qubits))]
    qubits = np.array(qubits)
    lent = int(len(qubits)/2)
    qubits.resize((lent,2))
    qubits = qubits.tolist()
    n = len(outputs)
    states = makeallstates(n)
    data = []
    for i in range(n):
        data.append(getdata(qubits,times,states[i]))
    fidels = [fidelity(data[k],outputs[k]) for k in range(n)]
    fidel = np.average(fidels)
    return -fidel

def solve(initqubits):
    res = minimize(calfil, initqubits, method='BFGS',options={'xtol': 1e-8, 'disp': True})
    return res






