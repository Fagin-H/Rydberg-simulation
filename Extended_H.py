# -*- coding: utf-8 -*-
"""
Created on Fri May 18 16:57:21 2018
@author: fhales
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt

I = qeye(2)
times = np.linspace(0, 7, 1000)
factors = np.linspace(1,0.1,100)
C3 = 7965
C6u = 766000
C6d = -5000

sigmap = sigmax() + sigmay()*1j
sigmam = sigmax() - sigmay()*1j

qubits = np.array([[0,0],[20,0],[40,0]])


def makematrix(qubits_co):
    intmatrix = [[np.linalg.norm(qubits_co[i]-qubits_co[j])**3 for j in range(len(qubits_co))] for i in range(len(qubits_co))]
    return intmatrix

def makeinputoutput(atom_number):
    tempstatein = fock(2,1)
    for i in range(atom_number-1):
        tempstatein = tensor(tempstatein,fock(2,0))
    tempstateout = tempstatein*tempstatein.dag()
    return tempstatein, tempstateout

def makesigf(i,j,atoms):
    temp1 = [I]*atoms
    temp2 = [I]*atoms
    temp1[i] = sigmap
    temp1[j] = sigmam
    temp2[i] = sigmam
    temp2[j] = sigmap
    temp = tensor(temp1)+tensor(temp2)
    return temp

def makesigz(i,j,ud,atoms):
    temp = [I]*atoms
    temp[i] = fock(2,ud)*fock(2,ud).dag()
    temp[j] = fock(2,ud)*fock(2,ud).dag()
    temp = tensor(temp)
    return temp

def makeham1(intmatrix):
    atoms = len(intmatrix)
    components =  []
    for i in range(atoms):
        for j in range(atoms):
            if j < i:
                components.append((C3*0.5/intmatrix[i][j]*makesigf(i,j,atoms)))
    ham = sum(components)
    return ham

def makeham2(intmatrix):
    atoms = len(intmatrix)
    components =  []
    for i in range(atoms):
        for j in range(atoms):
            if j < i:
                components.append((C3*0.5/intmatrix[i][j]*makesigf(i,j,atoms))+C6d/(intmatrix[i][j]**2)*makesigz(i,j,0,atoms)+C6u/(intmatrix[i][j]**2)*makesigz(i,j,1,atoms))
    ham = sum(components)
    return ham

def makeplot(H,timesteps,instate,basis):
    data = mcsolve(H,instate,timesteps,[],[])
    newdata=[]
    for i in range(len(timesteps)):
        newdata.append((fock(2,0).dag()*data.states[i].ptrace([len(instate.dims[0])-1])*fock(2,0)).tr())
    plt.plot(timesteps, newdata)

def doall(qubits,times):
    int_matrix = makematrix(qubits)
    input_state, output_basis = makeinputoutput(len(qubits))
    H = makeham1(int_matrix)
    makeplot(H,times,input_state,output_basis)
    
def make_entangled_inputoutput(atom_number):
    #For this we require, something (globally) like: qubits = np.array([[-100,0],[0,0],[20,0],[40,0]])
    tempstatein = bell_state(state='00')
    for i in range(atom_number-2):
        tempstatein = tensor(tempstatein,fock(2,0))
    tempstateout = tempstatein*tempstatein.dag()
    return tempstatein, tempstateout   

def makeplot_entanglement(H,timesteps,instate,basis):
    data = mcsolve(H,instate,timesteps,[],[])
    entanglement_concurrence = []
    for i in range(len(timesteps)):
       entanglement_concurrence.append(concurrence(data.states[i].ptrace([0,len(instate.dims[0])-1])))    
    plt.plot(timesteps, entanglement_concurrence)
    

def comph(qubits,times):
    input_state, output_basis = makeinputoutput(len(qubits))
    dists = []
    for i in range(len(factors)):
        int_matrix = makematrix(factors[i]*qubits)
        H1 = makeham1(int_matrix)
        H2 = makeham2(int_matrix)
        dists.append(comp1(H1,H2,input_state,(factors[i]**3)*times,output_basis))
    plt.plot(factors*20,dists)
    plt.ylabel('Graph difference')
    plt.xlabel('Atom separation (μm)')
    plt.savefig('Graph.png')
    
    
def comp1(H1,H2,input_state,timesx,output_basis):
    data1 = mcsolve(H1,input_state,timesx,[],output_basis)
    data2= mcsolve(H2,input_state,timesx,[],output_basis)
    data1.expect[0]
    return np.linalg.norm(data1.expect[0]-data2.expect[0])
    
    





















