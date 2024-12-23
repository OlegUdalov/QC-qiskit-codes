# Function used in VQE algorithms


import qiskit
from qiskit import QuantumRegister as Q_R
from qiskit import ClassicalRegister as C_R
import numpy as np

PI = np.pi

# Particles conserving anzatz
#This block of functions is used to create an Anzatz for VQE algorithm
#The anzatz can be used for modelling of quantum systems
#The anzatz conserves number of particles
#The anzatz includes single and double excitations
#The anzatz is taken from the paper "Efficient quatum circuits for quantum computational chemistry" by Yordanov
#See details on the anzatz in the  Jupyter notebook Cluster_anzatz

def u_se_1(theta, i, k, qc):
    '''
    This function creates a portion of single exciation operator
    Here we assume that i < k. 
    However, transitions goes both ways.
    '''
    PI = np.pi
    gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [k])
    qc.h(i)
    for i_qubit in range(k, i, -1):
        qc.cx(i_qubit, i_qubit - 1)
    gate = qiskit.circuit.library.RZGate(theta)
    qc.append(gate, [i])
    for i_qubit in range(i, k):
        qc.cx(i_qubit + 1, i_qubit)
    qc.h(i)
    gate = qiskit.circuit.library.RXGate(- PI / 2)
    qc.append(gate, [k])
    return qc

def u_se_2(theta, i, k, qc):
    '''
    The function realizes the circuit for the single exciattion term conjugated to the one above
    Here we assume that i < k for
    '''
    PI = np.pi
    gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [i])
    qc.h(k)
    for i_qubit in range(k, i, -1):
        qc.cx(i_qubit, i_qubit - 1)
    gate = qiskit.circuit.library.RZGate(-theta)
    qc.append(gate, [i])
    for i_qubit in range(i, k):
        qc.cx(i_qubit + 1, i_qubit)
    qc.h(k)
    gate = qiskit.circuit.library.RXGate(- PI / 2)
    qc.append(gate, [i])
    return qc

def u_se(theta, i, k, qc):
    '''
    Single excitation operator inducing transition between states i and k. Thetat defines how much of the state is transferring.
    Theta = 0 - no transition |i> -> |i>
    Theta = pi/2 - full transition  |i> -> |k>
    Theta = pi/4 - |i> -> 1/sqrt(2)(|i> + |k>)

    Parameters
    float Theta- angle of rotation. defines the strength of the transition
    int i, k - intial and final states
    QuantumCircuit qc - quantum circuit
    
    '''
    qc = u_se_1(theta, i, k, qc)
    qc = u_se_2(theta, i, k, qc)
    return qc

def staircase(theta, i, j, k, l, qc):
    '''
    Auxiliary function creating a staircase of cx gates with rotation.
    It is used in the double excitation operator
    i<j<k<l

    Parameters
    float Theta- angle of rotation gate in the middle of the staircase
    int i, j, k, l - intial and final states
    QuantumCircuit qc - quantum circuit
    
    '''
    PI = np.pi
    for i_q in range(l,k,-1):
        qc.cx(i_q, i_q - 1)
    qc.cx(k, j)
    for i_q in range(j,i,-1):
        qc.cx(i_q, i_q - 1)
    gate = qiskit.circuit.library.RZGate(- theta)
    qc.append(gate, [i])
    for i_q in range(i,j):
        qc.cx(i_q + 1, i_q)
    qc.cx(k, j)    
    for i_q in range(k,l):
        qc.cx(i_q + 1, i_q)
    return qc

def u_de(theta, i, j, k, l, qc):
    '''
    Creates an operator perfroming double excitation (transition from state |ij> to state |kl>)
    Theta = 0 - no transition |ij> -> |ij>
    Theta = pi/4 - full transition  |ij> -> |kl>
    Theta = pi/8 - |ij> -> 1/sqrt(2)(|ij> + |kl>)
    i<j, k<l

    Parameters
    float Theta- angle of rotation gate in the middle of the staircase
    int i, j, k, l - intial and final states
    QuantumCircuit qc - quantum circuit
    
    '''
    PI = np.pi
    if j < k:
        i_st = i
        j_st = j
        k_st = k
        l_st = l

    if j > k and j < l:
        i_st = i
        j_st = k
        k_st = j
        l_st = l

    if j > k and j > l:
        i_st = i
        j_st = k
        k_st = l
        l_st = j
    
    #1
    qc.h(k)
    qc.h(l)
    qc.h(i)
    gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [j])
    qc = staircase(theta, i_st, j_st, k_st, l_st, qc)
    qc.h(k)
    qc.h(l)
    qc.h(i)
    gate = qiskit.circuit.library.RXGate(- PI / 2)
    qc.append(gate, [j])

    #2
    qc.h(k)
    qc.h(l)
    qc.h(j)
    gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [i])
    qc = staircase(theta, i_st, j_st, k_st, l_st, qc)
    qc.h(k)
    qc.h(l)
    qc.h(j)
    gate = qiskit.circuit.library.RXGate(- PI / 2)
    qc.append(gate, [i])

    #3
    qc.h(l)
    gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [k])
    #gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [i])
    #gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [j])
    qc = staircase(theta, i_st, j_st, k_st, l_st, qc)
    qc.h(l)
    gate = qiskit.circuit.library.RXGate( - PI / 2)
    qc.append(gate, [k])
    #gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [i])
    #gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [j])

    #4
    qc.h(k)
    gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [l])
    #gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [i])
    #gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [j])
    qc = staircase(theta, i_st, j_st, k_st, l_st, qc)
    qc.h(k)
    gate = qiskit.circuit.library.RXGate( - PI / 2)
    qc.append(gate, [l])
    #gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [i])
    #gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [j])
    
    #5
    qc.h(l)
    gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [k])
    qc.h(i)
    qc.h(j)
    qc = staircase(-theta, i_st, j_st, k_st, l_st, qc)
    qc.h(l)
    gate = qiskit.circuit.library.RXGate(-PI / 2)
    qc.append(gate, [k])
    qc.h(i)
    qc.h(j)

    #6
    qc.h(k)
    gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [l])
    qc.h(i)
    qc.h(j)
    qc = staircase(-theta, i_st, j_st, k_st, l_st, qc)
    qc.h(k)
    gate = qiskit.circuit.library.RXGate(-PI / 2)
    qc.append(gate, [l])
    qc.h(i)
    qc.h(j)


    
    #7
    gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [l])
    qc.append(gate, [k])
    qc.h(j)
    qc.append(gate, [i])
    qc = staircase(-theta, i_st, j_st, k_st, l_st, qc)
    gate = qiskit.circuit.library.RXGate(-PI / 2)
    qc.append(gate, [l])
    qc.append(gate, [k])
    qc.h(j)
    qc.append(gate, [i])

    #8
    gate = qiskit.circuit.library.RXGate(PI / 2)
    qc.append(gate, [l])
    qc.append(gate, [k])
    qc.h(i)
    qc.append(gate, [j])
    qc = staircase(-theta, i_st, j_st, k_st, l_st, qc)
    gate = qiskit.circuit.library.RXGate(-PI / 2)
    qc.append(gate, [l])
    qc.append(gate, [k])
    qc.h(i)
    qc.append(gate, [j])
    
    return qc

def exc(n_qubits, theta, qc):
    '''
    The function creates an operator that perform single-particle and double-particles transitions on an initial state.
    Initial state should be defined outside of the function.
    The function is independent of the initial state.
    All pairs of single transitions |i> -> |j> are taken into account
    All pairs of double transitions |ij> -> |kl> are taken into account

    Parameters:
    int n_qubits - number of qubits in the cirquit qc. Should be more than 4. Corresponds to the number of states in the system.
    float array theta - array of angles theta, size of theta is defined as follows:
    len_se = n_qubits * (n_qubits - 1) / 2
    len_de = n_qubits * (n_qubits - 1) * (n_qubits - 2) * (n_qubits - 3) / 4 / 2
    len_tot = int(len_se + len_de)

    QuantumCircuit qc - quantum circuit
    
    '''
    n = 0
    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            qc = u_se(theta[n], i, j, qc)
            n = n + 1
    
    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            for k in range(i + 1, n_qubits):
                for l in range(k + 1, n_qubits):
                    if k != j and l !=j:
                        qc = u_de(theta[n], i, j, k, l, qc)
                        n = n + 1
    return qc
    
            