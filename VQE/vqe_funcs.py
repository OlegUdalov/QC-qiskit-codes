# Function used in VQE algorithms


import qiskit

from qiskit import QuantumRegister as Q_R
from qiskit import ClassicalRegister as C_R
from qiskit.quantum_info import SparsePauliOp

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


def ry_c(theta, i, k, qc):
    '''
    this is auxiliary function for single fermionic excitation. 
    This function is used to create Ansatz described in the paper 
    PHYSICAL REVIEW A 102, 062612 (2020)
    Efficient quantum circuits for quantum computational chemistry
    Yordan S. Yordanov, David R. M. Arvidsson-Shukur, and Crispin H. W. Barnes

    I call it cluster Ansatz by Yordanov. See corresponding Jupyter notebook about this ansatz

    Parameters:
    theta, float - rotation angle, defines how much of state i is transformed to k
    i,k, int - states between which transoition happens
    qc, QuantumCircuit - quantum circuit to which the ry_c set of gates will be added
    '''
    qc.cx(k,i)
    gate = qiskit.circuit.library.RZGate(PI / 2)
    qc.append(gate, [k])
    gate = qiskit.circuit.library.RYGate(-PI / 2)
    qc.append(gate, [i])
    gate = qiskit.circuit.library.RZGate(-PI / 2)
    qc.append(gate, [i])
    qc.cx(k,i)
    gate = qiskit.circuit.library.RYGate(theta / 2)
    qc.append(gate, [k])
    gate = qiskit.circuit.library.RZGate(-PI / 2)
    qc.append(gate, [i])
    qc.cx(k,i)
    gate = qiskit.circuit.library.RYGate(-theta / 2)
    qc.append(gate, [k])
    qc.h(i)
    #qc.cx(k,i)

    return qc

def ry_double_1(theta, i, j, k, l, qc):
    '''
    this is auxiliary function for double fermionic excitation

    This function is used to create Ansatz described in the paper 
    PHYSICAL REVIEW A 102, 062612 (2020)
    Efficient quantum circuits for quantum computational chemistry
    Yordan S. Yordanov, David R. M. Arvidsson-Shukur, and Crispin H. W. Barnes

    I call it cluster Ansatz by Yordanov. See corresponding Jupyter notebook about this ansatz

    Parameters:
    theta, float - rotation angle, defines how much of state i,j are transformed to k,l
    i,j,k,l int - states between which transoition happens
    qc, QuantumCircuit - quantum circuit to which the ry_c set of gates will be added
    '''
    #qc.cx(l, k)
    #qc.cx(j, i)
    qc.x(k)
    qc.x(i)
    #qc.cx(l, j)
    gate = qiskit.circuit.library.RYGate(theta / 8)
    qc.append(gate, [l])
    qc.h(k)
    qc.cx(l, k)
    gate = qiskit.circuit.library.RYGate(-theta / 8)
    qc.append(gate, [l])
    qc.h(i)
    qc.cx(l, i)
    gate = qiskit.circuit.library.RYGate(theta / 8)
    qc.append(gate, [l])
    qc.cx(l, k)
    gate = qiskit.circuit.library.RYGate(-theta / 8)
    qc.append(gate, [l])
    qc.h(j)
    qc.cx(l, j)
    gate = qiskit.circuit.library.RYGate(theta / 8)
    qc.append(gate, [l])
    qc.cx(l, k)
    gate = qiskit.circuit.library.RYGate(-theta / 8)
    qc.append(gate, [l])
    qc.cx(l, i)
    gate = qiskit.circuit.library.RYGate(theta / 8)
    qc.append(gate, [l])
    qc.h(i)
    qc.cx(l, k)
    gate = qiskit.circuit.library.RYGate(-theta / 8)
    qc.append(gate, [l])
    qc.h(k)
    qc.cx(l, j)
    qc.h(j)
    #gate = qiskit.circuit.library.RZGate(-PI / 2)
    #qc.append(gate, [j])
    #qc.cx(l, j)
    #gate = qiskit.circuit.library.RZGate(PI / 2)
    #qc.append(gate, [l])
    #gate = qiskit.circuit.library.RZGate(-PI / 2)
    #qc.append(gate, [j])
    qc.x(k)
    #gate = qiskit.circuit.library.RYGate(-PI / 2)
    #qc.append(gate, [j])
    qc.x(i)
    #qc.cx(l, k)
    #qc.cx(j, i)
    
    return qc

def se_yordanov(theta, i, k, qc):
    '''
    this is a function for single fermionic excitation
    Ansatz is described in the paper 
    PHYSICAL REVIEW A 102, 062612 (2020)
    Efficient quantum circuits for quantum computational chemistry
    Yordan S. Yordanov, David R. M. Arvidsson-Shukur, and Crispin H. W. Barnes

    I call it cluster Ansatz by Yordanov. See corresponding Jupyter notebook about this ansatz

    Parameters:
    theta, float - rotation angle, defines how much of state i is transformed to k, full transition is when theta = pi
    i,k, int - states between which transoition happens
    qc, QuantumCircuit - quantum circuit to which the ry_c set of gates will be added
    '''
    if k - 1 > i + 1: 
        for i_q in range(k - 1, i + 1, -1):
            qc.cx(i_q, i_q - 1)
    qc.cx(k, i)
    if  k - 1 > i + 1: 
        qc.cz(i + 1, k)
    ry_c(theta, i, k, qc)
    if k - 1 > i + 1: 
        qc.cz(i + 1, k)
    qc.cx(k, i)
    if k - 1 > i + 1: 
        for i_q in range(i + 1, k - 1):
            qc.cx(i_q + 1, i_q)
    
    return qc

def se_yordanov_no_ladder(theta, i, k, qc):
    '''
    this is a function for single fermionic excitation but without the ladder part
    Ansatz is described in the paper 
    PHYSICAL REVIEW A 102, 062612 (2020)
    Efficient quantum circuits for quantum computational chemistry
    Yordan S. Yordanov, David R. M. Arvidsson-Shukur, and Crispin H. W. Barnes

    I call it cluster Ansatz by Yordanov. See corresponding Jupyter notebook about this ansatz

    Parameters:
    theta, float - rotation angle, defines how much of state i is transformed to k, full transition is when theta = pi
    i,k, int - states between which transoition happens
    qc, QuantumCircuit - quantum circuit to which the ry_c set of gates will be added
    '''
    #if k - 1 > i + 1: 
    #    for i_q in range(k - 1, i + 1, -1):
    #        qc.cx(i_q, i_q - 1)
    qc.cx(k, i)
    if  k - 1 > i + 1: 
        qc.cz(i + 1, k)
    ry_c(theta, i, k, qc)
    if k - 1 > i + 1: 
        qc.cz(i + 1, k)
    qc.cx(k, i)
    #if k - 1 > i + 1: 
    #    for i_q in range(i + 1, k - 1):
    #        qc.cx(i_q + 1, i_q)
    
    return qc

def de_yordanov(theta, i, j, k, l, qc):
    '''
    this is a function for double fermionic excitation

    Ansatz is described in the paper 
    PHYSICAL REVIEW A 102, 062612 (2020)
    Efficient quantum circuits for quantum computational chemistry
    Yordan S. Yordanov, David R. M. Arvidsson-Shukur, and Crispin H. W. Barnes

    I call it cluster Ansatz by Yordanov. See corresponding Jupyter notebook about this ansatz

    Parameters:
    theta, float - rotation angle, defines how much of state i,j is transformed to k,l. Full transition is when theta = pi
    i,j,k,l, int - states between which transoition happens
    qc, QuantumCircuit - quantum circuit to which the ry_c set of gates will be added
    '''
    #theta = 8 * theta
    qc.cx(l, k)
    qc.cx(j, i)
    if l - 1 > i + 1: 
        for i_q in range(l - 1, i + 1, -1):
            qc.cx(i_q, i_q - 1)
    qc.cx(l, j)
    qc.cz(i + 1, l)
    
    ry_double_1(theta, i, j, k, l, qc)
    qc.cz(i + 1, l)
    qc.cx(l, j)
    if l - 1 > i + 1: 
        for i_q in range(i + 1, l - 1):
            qc.cx(i_q + 1, i_q)
    qc.cx(l, k)
    qc.cx(j, i)
    
    return qc

def de_yordanov_no_stair(theta, i, j, k, l, qc):
    '''
    this is a function for double qubit excitation. There is no staircase parts. So, the pajrity is not conerved and the operator does not make the right sign of the wave function. Absence of the staircase essentially reduces number of cnots in the circuit.

    Ansatz isea is taken from the paper 
    PHYSICAL REVIEW A 102, 062612 (2020)
    Efficient quantum circuits for quantum computational chemistry
    Yordan S. Yordanov, David R. M. Arvidsson-Shukur, and Crispin H. W. Barnes

    I call it cluster Ansatz by Yordanov without the staircase. See corresponding Jupyter notebook about this ansatz

    Parameters:
    theta, float - rotation angle, defines how much of state i,j is transformed to k,l. Full transition is when theta = pi
    i,j,k,l, int - states between which transoition happens
    qc, QuantumCircuit - quantum circuit to which the ry_c set of gates will be added
    '''
    #theta = 8 * theta
    qc.cx(l, k)
    qc.cx(j, i)
    qc.cx(l, j)
    qc.cz(i + 1, l)
    ry_double_1(theta, i, j, k, l, qc)
    qc.cz(i + 1, l)
    qc.cx(l, j)
    qc.cx(l, k)
    qc.cx(j, i)
    
    return qc


def sigma_x(node, nodes_number):
    if node > nodes_number - 1:
        print('Error: node index (''node'') cannot be higher than the total number of nodes')
        return -1
    interaction_string_1 = ''
    #term one
    if node > 0:
        for i in range(0, node):
            interaction_string_1 = interaction_string_1 + 'I'
    interaction_string_1 = interaction_string_1 + 'X'
    for i in range(node + 1, node + nodes_number):
        interaction_string_1 = interaction_string_1 + 'Z'
    interaction_string_1 = interaction_string_1 + 'Y'
    for i in range(node + nodes_number + 1, nodes_number * 2):
        interaction_string_1 = interaction_string_1 + 'I'
    #term two
    interaction_string_2 = ''
    if node > 0:
        for i in range(0, node):
            interaction_string_2 = interaction_string_2 + 'I'
    interaction_string_2 = interaction_string_2 + 'Y'
    for i in range(node + 1, node + nodes_number):
        interaction_string_2 = interaction_string_2 + 'Z'
    interaction_string_2 = interaction_string_2 + 'X'
    for i in range(node + nodes_number + 1, nodes_number * 2):
        interaction_string_2 = interaction_string_2 + 'I'
    bits = range(nodes_number * 2)
    interactions = [(interaction_string_1, bits, 1j/2)]
    interactions.append((interaction_string_2, bits, 1j/2))
    hamiltonian = SparsePauliOp.from_sparse_list(interactions, num_qubits = 2 * nodes_number)
    return hamiltonian

def sigma_y(node, nodes_number):
    if node > nodes_number - 1:
        print('Error: node index (''node'') cannot be higher than the total number of nodes')
        return -1
    interaction_string_1 = ''
    #term one
    if node > 0:
        for i in range(0, node):
            interaction_string_1 = interaction_string_1 + 'I'
    interaction_string_1 = interaction_string_1 + 'X'
    for i in range(node + 1, node + nodes_number):
        interaction_string_1 = interaction_string_1 + 'Z'
    interaction_string_1 = interaction_string_1 + 'X'
    for i in range(node + nodes_number + 1, nodes_number * 2):
        interaction_string_1 = interaction_string_1 + 'I'
    #term two
    interaction_string_2 = ''
    if node > 0:
        for i in range(0, node):
            interaction_string_2 = interaction_string_2 + 'I'
    interaction_string_2 = interaction_string_2 + 'Y'
    for i in range(node + 1, node + nodes_number):
        interaction_string_2 = interaction_string_2 + 'Z'
    interaction_string_2 = interaction_string_2 + 'Y'
    for i in range(node + nodes_number + 1, nodes_number * 2):
        interaction_string_2 = interaction_string_2 + 'I'
    bits = range(nodes_number * 2)
    interactions = [(interaction_string_1, bits, -1j/2)]
    interactions.append((interaction_string_2, bits, -1j/2))
    hamiltonian = SparsePauliOp.from_sparse_list(interactions, num_qubits = 2 * nodes_number)
    return hamiltonian

def sigma_z(node, nodes_number):
    if node > nodes_number - 1:
        print('Error: node index (''node'') cannot be higher than the total number of nodes')
        return -1

    #term two
    interaction_string_2 = ''
    if node > 0:
        for i in range(0, node):
            interaction_string_2 = interaction_string_2 + 'I'
    interaction_string_2 = interaction_string_2 + 'Z'
    for i in range(node + 1, 2 * nodes_number):
        interaction_string_2 = interaction_string_2 + 'I'

    #term four
    interaction_string_4 = ''
    for i in range(0, node + nodes_number):
        interaction_string_4 = interaction_string_4 + 'I'
    interaction_string_4 = interaction_string_4 + 'Z'
    for i in range(node + nodes_number + 1, 2 * nodes_number):
        interaction_string_4 = interaction_string_4 + 'I'

    
    bits = range(nodes_number * 2)
    interactions = [(interaction_string_2, bits, -1/2)]
    interactions.append((interaction_string_4, bits, 1/2))
    hamiltonian = SparsePauliOp.from_sparse_list(interactions, num_qubits = 2 * nodes_number)
    return hamiltonian

def Coulomb_on_site_ess(nodes_number, Uc):
    '''
    This function creates on-site Coulomb repulsion Hamiltonian. Only essential part of the hamiltonian n_up * n_down is present.
    The rest is just proportional to total particle numbers and does nothing for particle conserving anzats.

    Parameter:
    nodes_number, int - number of nodes in the system
    Uc, float - on-site Coulomb repulsion energy

    Output:
    hamiltonian
    
    '''
    interactions = []
    hamiltonian = []
    bits = range(nodes_number * 2)
    for i_node in range(0, nodes_number):
                   
        interaction_string_2 = ''
        if i_node > 0:
            for i in range(i_node):
                interaction_string_2 = interaction_string_2 + 'I'
        interaction_string_2 = interaction_string_2 + 'Z'
        for i in range(i_node + 1, i_node + nodes_number):
            interaction_string_2 = interaction_string_2 + 'I'
        interaction_string_2 = interaction_string_2 + 'Z'
        
        interactions.append((interaction_string_2, bits, Uc/2))
        
    hamiltonian =  (SparsePauliOp.from_sparse_list(interactions, num_qubits = 2 * nodes_number))
    return hamiltonian

def kinetic_energy(t, nodes_number, periodic = True):
    interactions = []
    hamiltonian = []
    bits = range(nodes_number * 2)
    for i_node in range(0, nodes_number - 1):
        interaction_string_1 = ''
        for i in range(i_node):
            interaction_string_1 = interaction_string_1 + 'I'
        interaction_string_1 = interaction_string_1 + 'XX'
        #print(interaction_string_1)
        interaction_string_2 = ''
        for i in range(i_node):
            interaction_string_2 = interaction_string_2 + 'I'
        interaction_string_2 = interaction_string_2 + 'YY'
        
        interaction_string_3 = ''
        for i in range(i_node + nodes_number):
            interaction_string_3 = interaction_string_3 + 'I'
        interaction_string_3 = interaction_string_3 + 'XX'
        
        interaction_string_4 = ''
        for i in range(i_node + nodes_number):
            interaction_string_4 = interaction_string_4 + 'I'
        interaction_string_4 = interaction_string_4 + 'YY'
        interactions.append((interaction_string_1, bits, -t/2))
        interactions.append((interaction_string_2, bits, -t/2))
        interactions.append((interaction_string_3, bits, -t/2))
        interactions.append((interaction_string_4, bits, -t/2))

    if periodic == True:
        interaction_string_1 = ''
        interaction_string_1 = interaction_string_1 + 'X'
        for i in range(1, nodes_number - 1):
            interaction_string_1 = interaction_string_1 + 'I'
        interaction_string_1 = interaction_string_1 + 'X'
        
        interaction_string_2 = ''
        interaction_string_2 = interaction_string_2 + 'Y'
        for i in range(1, nodes_number - 1):
            interaction_string_2 = interaction_string_2 + 'I'
        interaction_string_2 = interaction_string_2 + 'Y'
        
        interaction_string_3 = ''
        for i in range(0, nodes_number):
            interaction_string_3 = interaction_string_3 + 'I'
        interaction_string_3 = interaction_string_3 + 'X'
        for i in range(nodes_number + 1, 2 * nodes_number - 1):
            interaction_string_3 = interaction_string_3 + 'I'
        interaction_string_3 = interaction_string_3 + 'X'
        
        interaction_string_4 = ''
        for i in range(0, nodes_number):
            interaction_string_4 = interaction_string_4 + 'I'
        interaction_string_4 = interaction_string_4 + 'Y'
        for i in range(nodes_number + 1, 2 * nodes_number - 1):
            interaction_string_4 = interaction_string_4 + 'I'
        interaction_string_4 = interaction_string_4 + 'Y'
        
        interactions.append((interaction_string_1, bits, -t/2))
        interactions.append((interaction_string_2, bits, -t/2))
        interactions.append((interaction_string_3, bits, -t/2))
        interactions.append((interaction_string_4, bits, -t/2))

    hamiltonian =  (SparsePauliOp.from_sparse_list(interactions, num_qubits = 2 * nodes_number))
    return hamiltonian


def full_ham(nodes_number, magnetization, J, t, Uc, periodic = True):
    '''
    Create 1D Hamiltonian consisting of tight-binding kinetic energy (with energy constant t), 
    periodic boundary conditiona may be applied
    s-d interaction with local magnetic moments and interaction conatnt J

    Parameters:
    nodes_number, int - number of nodes in the system
    magnetization, 2D float array [nodes_number x 3] - magnetization direction for each node
    J, float - s-d exchange constant
    t, float - kinetic energy constant
    periodic, boolean - True means that periodic boundary conditions are applied
    
    '''
    hamiltonian = []
    hamiltonian = kinetic_energy(t, nodes_number, periodic)
    #mag_ham = SparsePauliOp.from_sparse_list([], num_qubits = 2 * nodes_number)
    #mag_ham = sum([], J * magnetization[1][2] * sigma_x(1, nodes_number))
    mag_ham = []
    for node in range(nodes_number):
        if magnetization[node][0] != 0:
            #mag_ham.append(J*magnetization[node][0] * sigma_x(node, nodes_number))
            mag_ham = sum(mag_ham, J * magnetization[node][0] * sigma_x(node, nodes_number))
        if magnetization[node][1] != 0:
            #mag_ham.append(J*magnetization[node][1] * sigma_y(node, nodes_number))
            mag_ham = sum(mag_ham, J * magnetization[node][1] * sigma_y(node, nodes_number))
        if magnetization[node][2] != 0:
            #mag_ham.append(J*magnetization[node][2] * sigma_z(node, nodes_number))
            mag_ham = sum(mag_ham, J * magnetization[node][2] * sigma_z(node, nodes_number))
    hamiltonian = sum(mag_ham, hamiltonian)
    coul_ham = Coulomb_on_site_ess(nodes_number, Uc)
    hamiltonian = sum(coul_ham, hamiltonian)
    return hamiltonian

def exc_yordanov(n_qubits, theta, qc):
    '''
    The function creates an operator that perform single-particle and double-particles transitions on an initial state.
    Initial state should be defined outside of the function.
    The function is independent of the initial state.
    
    Ansatz is described in the paper 
    PHYSICAL REVIEW A 102, 062612 (2020)
    Efficient quantum circuits for quantum computational chemistry
    Yordan S. Yordanov, David R. M. Arvidsson-Shukur, and Crispin H. W. Barnes

    I call it cluster Ansatz by Yordanov. See corresponding Jupyter notebook about this ansatz
    

    Parameters:
    n_qubits - number of qubits in the cirquit qc. Should be more than 4. Corresponds to the number of states in the system.
    theta - array of angles theta, 0<theta<pi 
    
    '''
    n = 0
    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            qc = se_yordanov(theta[n], i, j, qc)
            n = n + 1
    
    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            for k in range(i + 1, n_qubits):
                for l in range(k + 1, n_qubits):
                    if k != j and l !=j:
                        qc = de_yordanov(theta[n], i, j, k, l, qc)
                        n = n + 1
    return qc

def exc_yordanov_no_stair(n_qubits, theta, qc):
    '''
    The function creates an operator that perform single-particle and double-particles transitions on an initial state.
    Initial state should be defined outside of the function.
    The function is independent of the initial state.
    
    Ansatz is described in the paper 
    PHYSICAL REVIEW A 102, 062612 (2020)
    Efficient quantum circuits for quantum computational chemistry
    Yordan S. Yordanov, David R. M. Arvidsson-Shukur, and Crispin H. W. Barnes

    I call it cluster Ansatz by Yordanov. See corresponding Jupyter notebook about this ansatz
    

    Parameters:
    n_qubits - number of qubits in the cirquit qc. Should be more than 4. Corresponds to the number of states in the system.
    theta - array of angles theta, 0<theta<pi 
    
    '''
    n = 0
    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            qc = se_yordanov(theta[n], i, j, qc)
            n = n + 1
    
    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            for k in range(i + 1, n_qubits):
                for l in range(k + 1, n_qubits):
                    if k != j and l !=j:
                        qc = de_yordanov_no_stair(theta[n], i, j, k, l, qc)
                        n = n + 1
    return qc

def exc_yordanov_single_only(n_qubits, theta, qc, n_red):
    '''
    The function creates an operator that perform single-particle transitions only on an initial state.
    Initial state should be defined outside of the function.
    The function is independent of the initial state.
    
    Ansatz is described in the paper 
    PHYSICAL REVIEW A 102, 062612 (2020)
    Efficient quantum circuits for quantum computational chemistry
    Yordan S. Yordanov, David R. M. Arvidsson-Shukur, and Crispin H. W. Barnes

    I call it cluster Ansatz by Yordanov. See corresponding Jupyter notebook about this ansatz
    

    Parameters:
    n_qubits - number of qubits in the cirquit qc. Should be more than 4. Corresponds to the number of states in the system.
    theta - array of angles theta, 0<theta<pi 
    
    '''
    n = 0
    for i in range(n_qubits - n_red):
        for j in range(i + 1, n_qubits):
            qc = se_yordanov(theta[n], i, j, qc)
            n = n + 1
    '''
    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            for k in range(i + 1, n_qubits):
                for l in range(k + 1, n_qubits):
                    if k != j and l !=j:
                        qc = de_yordanov_no_stair(theta[n], i, j, k, l, qc)
                        n = n + 1
    '''
    return qc

def exc_yordanov_single_only_no_ladder(n_qubits, theta, qc, n_red):
    '''
    The function creates an operator that perform single-particle transitions only on an initial state.
    Initial state should be defined outside of the function.
    The function is independent of the initial state.
    
    Ansatz is described in the paper 
    PHYSICAL REVIEW A 102, 062612 (2020)
    Efficient quantum circuits for quantum computational chemistry
    Yordan S. Yordanov, David R. M. Arvidsson-Shukur, and Crispin H. W. Barnes

    I call it cluster Ansatz by Yordanov. See corresponding Jupyter notebook about this ansatz
    

    Parameters:
    n_qubits - number of qubits in the cirquit qc. Should be more than 4. Corresponds to the number of states in the system.
    theta - array of angles theta, 0<theta<pi 
    
    '''
    n = 0
    for i in range(n_qubits - n_red - 1):
        for j in range(i + 1, n_qubits):
            qc = se_yordanov_no_ladder(theta[n], i, j, qc)
            n = n + 1
    
    '''
    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            for k in range(i + 1, n_qubits):
                for l in range(k + 1, n_qubits):
                    if k != j and l !=j:
                        qc = de_yordanov_no_stair(theta[n], i, j, k, l, qc)
                        n = n + 1
    '''
    return qc