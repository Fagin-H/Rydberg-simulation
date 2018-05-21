# -*- coding: utf-8 -*-
"""
Created on Fri May 18 16:57:21 2018
@author: fhales
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt

I = qeye(2)
times = np.linspace(0.0, 100.0, 1000.0)
C3=1


qubits = [[0,0],[1,0],[2,0]]

def makematrix(qubits_co):
    intmatrix = []
    for i in range(len(qubits_co)):
        temprow = []
        for j in range(len(qubits_co)):
            temprow.append(((qubits_co[i][0]-qubits_co[j][0])**2+(qubits_co[i][1]-qubits_co[j][1])**2)**1.5)
        intmatrix.append(temprow)
    return intmatrix
        
def makeinputoutput(atom_number):
    tempstatein = fock(2,1)
    for i in range(atom_number-1):
        tempstatein = tensor(tempstatein,fock(2,0))
    tempstateout = tempstatein*tempstatein.dag()
    return tempstatein, tempstateout

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
    ham = 0.5*C3*sum(components)
    return ham


def makeplot(H,timesteps,instate,basis):
    data = mcsolve(H,instate,timesteps,[],basis)
    plt.plot(timesteps, data.expect[0])

def doall(qubits,times):
    int_matrix = makematrix(qubits)
    input_state, output_basis = makeinputoutput(len(qubits))
    H = makeham(int_matrix)
    makeplot(H,times,input_state,output_basis)