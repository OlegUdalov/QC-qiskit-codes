'''this file contains some quantum routines that I will use in Shor algorithm'''


import qiskit
import numpy as np
import math

from qiskit_aer import Aer
from qiskit_aer import AerSimulator
from qiskit import QuantumCircuit, transpile
from qiskit.visualization import plot_histogram
from qiskit.circuit.library import SwapGate
from qiskit.circuit.library import PhaseGate
from qiskit import QuantumRegister as Q_R
from qiskit import ClassicalRegister as C_R
import classical_routines as cr

PI=np.pi

def QFTn_meas(n, psi_init, measure):
    """the function implements QFT transformation QC with n qubits
    
    The algorithm follows the one provided at https://github.com/Qiskit/textbook/blob/main/notebooks/ch-algorithms/quantum-fourier-transform.ipynb 
    The function adds QFT circuit to the existing circuit circ using a list of qubits qreg
    
        
    Parameter:
    __________
    n: integer, number of qubits that will be used for QFT
    psi_init: list, should have 2^n numbers between 0 and 1, represents the state of n input qubits
    measure: boolean, if 'yes' the classical register is added to the circuit with n cells, measurements of the qubits are added after the QFT
    
    
    Returns:
    __________
    Depends on the measurement parameter
    if 'no'
    [QFT,qubits1]: a list, new quantum circuit with added QFT and the qubits1 register created inside of the function
    if 'yes'
    [QFT,qubits1,readout]: a list, new quantum circuit with added QFT, the qubits1 register created inside of the function, cleassical readout afte the 
    measurements

    Notes:
    __________
    Qiskit 1.0 was used to create the function

    References:
    ___________
    My primitive study of the QFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    qubits_number=n
    qubits1 = qiskit.QuantumRegister(qubits_number)
    if measure=='yes':
        Readout=qiskit.ClassicalRegister(qubits_number)
        QFT = qiskit.QuantumCircuit(qubits1,Readout)
    else:
        QFT = qiskit.QuantumCircuit(qubits1)
    
    norm=np.linalg.norm(psi_init)
    psi_init=psi_init/norm
    QFT.initialize(psi_init)

    for i in range(0,qubits_number):
        QFT.h(qubits1[qubits_number-1-i])
        if i<QubitsNumber-1:
            for j in range(0,qubits_number-i-1):
                Rangle=PI/(math.pow(2,j+1))
                QFT.cp(Rangle,qubits1[qubits_number-1-i],qubits1[qubits_number-2-j-i])
    for i in range(0,qubits_number//2):
        QFT.swap(qubits1[i],qubits1[qubits_number-1-i])
    
    if measure=='yes':
        QFT.measure(qubits1, readout)
        return [QFT,qubits1,readout]
    else:
        return [QFT,qubits1]

def QFTn(qreg, circ):
    """the function implements QFT transformation QC with n qubits
    The algorithm follows the one provided at https://github.com/Qiskit/textbook/blob/main/notebooks/ch-algorithms/quantum-fourier-transform.ipynb 
    The function adds QFT circuit to the existing circuit circ using a list of qubits qreg
        
    Parameter:
    __________
    qreg: qiskit.QuantumRegister, list of qubits that will be used for QFT
    circ: qiskit.QuantumCircuit, a quantum circuit to which the QFT is added
    
    
    Returns:
    __________
    qiskit.QuantumCircuit, a new quantum circuit with added QFT


    Notes:
    __________
    Qiskit 1.0 was used to create the function

    References:
    ___________
    My primitive study of the QFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    qubits_number = len(qreg)

    for i in range(0, qubits_number):
        circ.h(qreg[qubits_number - 1 - i])
        if i<qubits_number - 1:
            for j in range(0, qubits_number - i - 1):
                r_angle = PI / (math.pow(2, j + 1))
                circ.cp(r_angle, qreg[qubits_number - 1 - i], qreg[qubits_number - 2 - j - i])
    for i in range(0, qubits_number // 2):
        circ.swap(qreg[i], qreg[qubits_number - 1 - i])
    
    return circ

def QFTn_instr(bit_size):
    """the function implements QFT transformation QC with bit_szie qubits
    The algorithm follows the one provided at https://github.com/Qiskit/textbook/blob/main/notebooks/ch-algorithms/quantum-fourier-transform.ipynb 
    The function creates instruction that can be added to an existing circuit
    The instruction cen be transformed to a gate. This gate have bit_size input (ant output) qubits and no classical register
        
    Parameter:
    __________
    bit_size: int, number of qubits for QFT
    
    
    
    Returns:
    __________
    qiskit.instruction, a new quantum instruction with QFT


    Notes:
    __________
    Qiskit 1.0 was used to create the function

    References:
    ___________
    My primitive study of the QFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    qreg = qiskit.QuantumRegister(bit_size) #this quantum register stores input quantum data
    circ = qiskit.QuantumCircuit(qreg)

    for i in range(0, bit_size):
        circ.h(qreg[bit_size - 1 - i])
        if i < bit_size - 1:
            for j in range(0, bit_size - i - 1):
                r_angle = PI / (math.pow(2, j + 1))
                circ.cp(r_angle, qreg[bit_size - 1 - i], qreg[bit_size - 2 - j - i])
    for i in range(0, bit_size // 2):
        circ.swap(qreg[i], qreg[bit_size - 1 - i])
    
    instr = circ.to_instruction(None, 'QFT\n size '+ str(bit_size))
    return instr


def QFTn_contr_gate(bit_size, contr_qubits_number):
    """the function implements QFT transformation QC with n qubits controlled by contr_qubits_number qubits
    The QFT algorithm itself  follows the one provided at https://github.com/Qiskit/textbook/blob/main/notebooks/ch-algorithms/quantum-fourier-transform.ipynb 
    The function create a controlled gate

    The gate has bit_size + contr_qubits_number inputs and same number of outputs
    First contr_qubits_number qubits are for the control
    The last bit_size qubits are the qubits for the QFT 
    
        
    Parameter:
    __________
    bit_size: int, number of qbits for QFT
    contr_qubits_number: number of control qubits
    
    Returns:
    __________
    qiskit.QuantumCircuit.gate, a gate with controlled QFT, first contr_qubits_number


    Notes:
    __________
    Qiskit 1.0 was used to create the function

    """
    qreg = qiskit.QuantumRegister(bit_size)
    circ = qiskit.QuantumCircuit(qreg)
    
    for i in range(0, bit_size):
        circ.h(qreg[bit_size - 1 - i])
        if i<bit_size - 1:
            for j in range(0, bit_size - i - 1):
                r_angle = PI / (math.pow(2, j + 1))
                circ.cp(r_angle, qreg[bit_size - 1 - i], qreg[bit_size - 2 - j - i])
    for i in range(0, bit_size // 2):
        circ.swap(qreg[i], qreg[bit_size - 1 - i])

    if contr_qubits_number > 0:
        qft_gate = circ.to_gate(None,'QFT').control(contr_qubits_number)
    else:
        qft_gate = circ.to_gate(None,'QFT')
    
    return qft_gate

def IQFTn(qreg, circ):
    """the function implements inverse QFT transformation QC with n qubits
    The algorithm follows the one provided at https://github.com/Qiskit/textbook/blob/main/notebooks/ch-algorithms/quantum-fourier-transform.ipynb 
    The function adds QFT circuit to the existing circuit circ using a list of qubits qreg
        
    Parameter:
    __________
    bit_size: int, number of qbits for QFT
    contr_qubits_number: number of control qubits
    
    
    Returns:
    __________
    qiskit.QuantumCircuit.gate, a new quantum circuit with added QFT


    Notes:
    __________
    Qiskit 1.0 was used to create the function

    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    qubits_number = len(qreg)
    for i in range(0, qubits_number // 2):
        circ.swap(qreg[i], qreg[qubits_number - 1 - i])
    
    for i in range(0, qubits_number):
        if i > 0:
            for j in range(0, i):
                r_angle = -PI / (math.pow(2, i - j))
                circ.cp(r_angle, qreg[i], qreg[j])
        circ.h(qreg[i])
    return circ

def IQFTn_instr(bit_size):
    """the function implements inverse QFT transformation QC with bit_size qubits
    The algorithm follows the one provided at https://github.com/Qiskit/textbook/blob/main/notebooks/ch-algorithms/quantum-fourier-transform.ipynb 
    The function acreates instruction
    The instruction cen be transformed to a gate. This gate have bit_size input (ant output) qubits and no classical register
        
    Parameter:
    __________
    bit_size: int, number of qbits for IQFT
       
    
    Returns:
    __________
    qiskit.QuantumCircuit.instruction, a quantum instruction


    Notes:
    __________
    Qiskit 1.0 was used to create the function

    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    qreg = qiskit.QuantumRegister(bit_size) #this quantum register stores input quantum data
    circ = qiskit.QuantumCircuit(qreg)
    
    for i in range(0, bit_size // 2):
        circ.swap(qreg[i], qreg[bit_size - 1 - i])
    
    for i in range(0, bit_size):
        if i > 0:
            for j in range(0, i):
                r_angle = -PI / (math.pow(2, i - j))
                circ.cp(r_angle, qreg[i], qreg[j])
        circ.h(qreg[i])
    instr = circ.to_instruction(None, 'IQFT\n size '+ str(bit_size))

    return instr

def IQFTn_contr_gate(bit_size, contr_qubits_number):
    """the function implements inverse QFT transformation QC with bit_size qubits. This inverse QFT is controlled by contr_qubits_number qubits
    The algorithm follows the one provided at https://github.com/Qiskit/textbook/blob/main/notebooks/ch-algorithms/quantum-fourier-transform.ipynb 
    The creates a controlled gate
    
    The gate has bit_size + contr_qubits_number inputs and same number of outputs
    First contr_qubits_number qubits are for the control
    The last bit_size qubits are the qubits for the IQFT 
        
    Parameter:
    __________
    qreg: qiskit.QuantumRegister, list of qubits that will be used for QFT
    circ: qiskit.QuantumCircuit, a quantum circuit to which the QFT is added
    
    
    Returns:
    __________
    qiskit.QuantumCircuit, a new quantum circuit with added QFT


    Notes:
    __________
    Qiskit 1.0 was used to create the function

    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    qreg = qiskit.QuantumRegister(bit_size)
    circ = qiskit.QuantumCircuit(qreg)
    for i in range(0, bit_size // 2):
        circ.swap(qreg[i], qreg[bit_size - 1 - i])
    
    for i in range(0, bit_size):
        if i > 0:
            for j in range(0, i):
                r_angle = -PI / (math.pow(2, i - j))
                circ.cp(r_angle, qreg[i], qreg[j])
        circ.h(qreg[i])
    
    if contr_qubits_number > 0:
        iqft_gate = circ.to_gate(None,'IQFT').control(contr_qubits_number)
    else:
        iqft_gate = circ.to_gate(None,'IQFT')
    
    return iqft_gate

def draper_adder_cl(q_reg, cl_reg, circ):
    """the function implements the simplified Draper adder which adds a quantum number to a clasical number
    The algorithm follows the paper arXiv:quant-ph/0205095
    The function creates a quantum circuit for Drapper-kind adder and adds this quantum circuit to the end of the existing one transferred as a prarameter 
    to this function. The adder adds a classical number given as a parameter cl_reg to a quantum number given as a parameter q_reg. This is addition of two 
    binary numbers.

    If number of classical bits is not equal to number of qubits, the function does nothing and return message.

    Parameter:
    __________
    q_reg: qiskit.QuantumRegister, list of qubits that will be used for addition (the first qubit in the list represents the highest digit of quantum number)
    cl_reg: list, list of classical bits that will be used for addition (the first bit in the list represents the highest digit of classical number)
    circ: qiskit.QuantumCircuit, a quantum circuit to which the QFT is added
    
    
    Returns:
    __________
    qiskit.QuantumCircuit, a new quantum circuit with added Draper adder
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    
    bit_size = len(q_reg) #size (number of digits) of the binary number
    cl_reg_size = len(cl_reg)
    if bit_size != cl_reg_size:
        print('draper_adder_cl: Size of quantum register is not the same as that of classical one, so the function does nothing')
 

    for i in range(0, int(np.floor(bit_size / 2))):
        circ.swap(q_reg[i], q_reg[bit_size - 1 - i])
        
    
    circ=QFTn(q_reg, circ)
    
    for i_q in range(bit_size): #loop ober the quantum register
        for i_r in range(bit_size - i_q):
            with circ.if_test((cl_reg[i_r + i_q], 1)):
                circ.p(2 * PI / (pow(2, i_r + 1)), q_reg[i_q])
    
    circ = IQFTn(q_reg, circ)

    for i in range(0, int(np.floor(bit_size / 2))):
        circ.swap(q_reg[i], q_reg[bit_size - 1 - i])

    return circ

def draper_adder_cl_instr(bit_size, control_bit_size):
    """the function implements the simplified Draper adder which adds a quantum number to a clasical number, the adder is controlled
    The algorithm follows the paper arXiv:quant-ph/0205095
    
    The adder is controlled with a number (control_bit_size) of qubits
        There are control_bit_size qubits performing control of the Draper adder. They both should have |1> state to turn on the adder.
    Control qubits are the upper control_bit_size quibits, added quantum number is the qubits below the control qubits
    
    The function creates an instruction (single block) having bit_size + control_bit_size qubits and bit_size classical bits as input and as output

    Parameter:
    __________
    bit_size, size of the binary numbers we add with this adder
    control_bit_size, number of control qubits
    
    Returns:
    __________
    qiskit.QuantumCircuit.instruction, a new instruction with quantumly controlled Draper adder
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    References:
    ___________
    My primitive study of the Draper adder can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    
    q_reg = qiskit.QuantumRegister(bit_size + control_bit_size) #this quantum register stores the added quantum number a
    cl_reg = qiskit.ClassicalRegister(bit_size) #this classical register stores classical number b
    adder = qiskit.QuantumCircuit(q_reg, cl_reg)

    if control_bit_size > 0:
        for i in range(0,int(np.floor(bit_size / 2))):
            aux_gate = SwapGate().control(control_bit_size)
            qbit_list = list(range(0, control_bit_size))
            qbit_list.append(i + control_bit_size)
            qbit_list.append(bit_size + control_bit_size - 1 - i)
            adder.append(aux_gate,qbit_list)
    else:
        for i in range(0,int(np.floor(bit_size / 2))):
            aux_gate = SwapGate()
            qbit_list = [i + control_bit_size]
            qbit_list.append(bit_size + control_bit_size - 1 - i)
            adder.append(aux_gate,qbit_list)        
    
   
    aux_gate = QFTn_contr_gate(bit_size, control_bit_size)
    adder.append(aux_gate,list(range(0, bit_size + control_bit_size)))
    
    if control_bit_size > 0:
        for i_q in range(bit_size): #loop over the quantum register
            for i_r in range(bit_size - i_q):
                aux_gate = PhaseGate(2 * PI / (pow(2, i_r + 1))).control(control_bit_size).c_if(cl_reg[i_r + i_q], 1)
                qbit_list = list(range(0,control_bit_size))
                qbit_list.append(i_q + control_bit_size)
                adder.append(aux_gate, qbit_list)
    else:
        for i_q in range(bit_size): #loop over the quantum register
            for i_r in range(bit_size - i_q):
                aux_gate = PhaseGate(2 * PI / (pow(2, i_r + 1))).c_if(cl_reg[i_r + i_q], 1)
                qbit_list = [i_q + control_bit_size]
                adder.append(aux_gate, qbit_list)        
    
    aux_gate = IQFTn_contr_gate(bit_size, control_bit_size)
    adder.append(aux_gate,list(range(bit_size + control_bit_size)))

    if control_bit_size > 0:
        for i in range(0,int(np.floor(bit_size / 2))):
            aux_gate = SwapGate().control(control_bit_size)
            qbit_list = list(range(0, control_bit_size))
            qbit_list.append(i + control_bit_size)
            qbit_list.append(bit_size + control_bit_size - 1 - i)
            adder.append(aux_gate,qbit_list)
    else:
        for i in range(0,int(np.floor(bit_size / 2))):
            aux_gate = SwapGate()
            qbit_list = [i + control_bit_size]
            qbit_list.append(bit_size + control_bit_size - 1 - i)
            adder.append(aux_gate,qbit_list)        
    
    instr = adder.to_instruction(None, 'ctrl_drap_add_cl\n ctr '+ str(control_bit_size))

    return instr

def phi_add_instr(bit_size, control_bit_size):
    """the function implements the simplified Draper adder which adds a QFT of quantum number to a clasical number, the adder is controlled
    The input should QFT of the quantum number a, the output is QFT(a+b), 
    The algorithm follows the paper arXiv:quant-ph/0205095
    
    The adder is controlled with a number (control_bit_size) of qubits
    There are control_bit_size qubits performing control of the Draper adder. They both should have |1> state to turn on the adder.
    Control qubits are the upper control_bit_size quibits, added quantum number is the qubits below the control qubits
    
    The function creates an instruction (single block) having bit_size + control_bit_size qubits and bit_size classical bits as input and as output

    Parameter:
    __________
    bit_size, size of the binary numbers we add with this adder
    control_bit_size, number of control qubits
    
    Returns:
    __________
    qiskit.QuantumCircuit.instruction, a new instruction with quantumly controlled Draper adder
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    References:
    ___________
    My primitive study of the Draper adder can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    
    q_reg = qiskit.QuantumRegister(bit_size + control_bit_size) #this quantum register stores the added quantum number a
    cl_reg = qiskit.ClassicalRegister(bit_size) #this classical register stores classical number b
    adder = qiskit.QuantumCircuit(q_reg, cl_reg)

    if control_bit_size > 0:
        for i_q in range(bit_size): #loop over the quantum register
            for i_r in range(bit_size - i_q):
                aux_gate = PhaseGate(2 * PI / (pow(2, i_r + 1))).control(control_bit_size).c_if(cl_reg[i_r + i_q], 1)
                qbit_list = list(range(0, control_bit_size))
                qbit_list.append(i_q + control_bit_size)
                adder.append(aux_gate, qbit_list)
    else:
        for i_q in range(bit_size): #loop over the quantum register
            for i_r in range(bit_size - i_q):
                aux_gate = PhaseGate(2 * PI / (pow(2, i_r + 1))).c_if(cl_reg[i_r + i_q], 1)
                qbit_list = [i_q + control_bit_size]
                adder.append(aux_gate, qbit_list)       
    
   
    instr = adder.to_instruction(None, 'ctrl_phi_add\n ctr '+ str(control_bit_size))

    return instr

def draper_subtraction_cl_instr(bit_size, control_bit_size):
    """the function implements the simplified Draper subtractor which subtracts a clasical number from quantum number 
    The algorithm follows the paper arXiv:quant-ph/0205095
    The function creates an instruction (single block) having n_bits qubits and n_bits classical bits as input and as output

    Parameter:
    __________
    bit_size, size of the binary numbers we add with this adder
    
    
    Returns:
    __________
    qiskit.QuantumCircuit,instruction, a new instruction with added Draper adder
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    
    
    q_reg = qiskit.QuantumRegister(bit_size + control_bit_size) #this quantum register stores the added quantum number a
    cl_reg = qiskit.ClassicalRegister(bit_size) #this classical register stores classical number b
    adder = qiskit.QuantumCircuit(q_reg, cl_reg)

    if control_bit_size > 0:
        for i in range(0,int(np.floor(bit_size / 2))):
            aux_gate = SwapGate().control(control_bit_size)
            qbit_list = list(range(0, control_bit_size))
            qbit_list.append(i + control_bit_size)
            qbit_list.append(bit_size + control_bit_size - 1 - i)
            adder.append(aux_gate,qbit_list)
    else:
        for i in range(0,int(np.floor(bit_size / 2))):
            aux_gate = SwapGate()
            qbit_list = [i + control_bit_size]
            qbit_list.append(bit_size + control_bit_size - 1 - i)
            adder.append(aux_gate,qbit_list)        

    aux_gate = IQFTn_contr_gate(bit_size, control_bit_size)
    adder.append(aux_gate,list(range(0, bit_size + control_bit_size)))
    
    if control_bit_size > 0:
        for i_q in range(bit_size): #loop over the quantum register
            for i_r in range(bit_size - i_q):
                aux_gate = PhaseGate(2 * PI / (pow(2, i_r + 1))).control(control_bit_size).c_if(cl_reg[i_r + i_q], 1)
                qbit_list = list(range(0,control_bit_size))
                qbit_list.append(i_q + control_bit_size)
                adder.append(aux_gate, qbit_list)
    else:
        for i_q in range(bit_size): #loop over the quantum register
            for i_r in range(bit_size - i_q):
                aux_gate = PhaseGate(2 * PI / (pow(2, i_r + 1))).c_if(cl_reg[i_r + i_q], 1)
                qbit_list = [i_q + control_bit_size]
                adder.append(aux_gate, qbit_list)
    
    aux_gate = QFTn_contr_gate(bit_size, control_bit_size)
    adder.append(aux_gate,list(range(bit_size + control_bit_size)))

    if control_bit_size > 0:
        for i in range(0,int(np.floor(bit_size / 2))):
            aux_gate = SwapGate().control(control_bit_size)
            qbit_list = list(range(0, control_bit_size))
            qbit_list.append(i + control_bit_size)
            qbit_list.append(bit_size + control_bit_size - 1 - i)
            adder.append(aux_gate,qbit_list)
    else:
        for i in range(0,int(np.floor(bit_size / 2))):
            aux_gate = SwapGate()
            qbit_list = [i + control_bit_size]
            qbit_list.append(bit_size + control_bit_size - 1 - i)
            adder.append(aux_gate,qbit_list)   
    
    adder.inverse()
    
    instr = adder.to_instruction(None, 'ctrl_draper_subtr_cl \n ctrl' + str(control_bit_size))

    return instr

def phi_subtr_instr(bit_size, control_bit_size):
    """the function implements the simplified Draper subtractor which subtracts a clasical number a from of quantum number b
    The input should QFT of the quantum number a, the output is QFT(a-b), 
    The algorithm follows the paper arXiv:quant-ph/0205095
    The function creates an instruction (single block) having n_bits qubits and n_bits classical bits as input and as output

    Parameter:
    __________
    bit_size, size of the binary numbers we add with this adder
    
    
    Returns:
    __________
    qiskit.QuantumCircuit,instruction, a new instruction with added Draper adder
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    
    
    q_reg = qiskit.QuantumRegister(bit_size + control_bit_size) #this quantum register stores the added quantum number a
    cl_reg = qiskit.ClassicalRegister(bit_size) #this classical register stores classical number b
    adder = qiskit.QuantumCircuit(q_reg, cl_reg)

    if control_bit_size > 0:
        for i_q in range(bit_size): #loop over the quantum register
            for i_r in range(bit_size - i_q):
                aux_gate = PhaseGate(2 * PI / (pow(2, i_r + 1))).control(control_bit_size).c_if(cl_reg[i_r + i_q], 1)
                qbit_list = list(range(0,control_bit_size))
                qbit_list.append(i_q + control_bit_size)
                adder.append(aux_gate, qbit_list)
    else:
        for i_q in range(bit_size): #loop over the quantum register
            for i_r in range(bit_size - i_q):
                aux_gate = PhaseGate(2 * PI / (pow(2, i_r + 1))).c_if(cl_reg[i_r + i_q], 1)
                qbit_list = [i_q + control_bit_size]
                adder.append(aux_gate, qbit_list)
    
    adder.inverse()
    
    instr = adder.to_instruction(None, 'phi_subtr \n ctrl' + str(control_bit_size))

    return instr

def ctrl_add_mod_N(bit_size, ctrl_bit_size):
    """the function implements the adder modulo N (N is classical number) which adds a clasical number a and quantum number b
    The algorithm follows the paper arXiv:quant-ph/0205095
    The function creates an instruction (single block) having bit_size + control_bit_size qubits and 2 * bit_size classical bits as input and as output
    Each number has bit_size bits
    The adder is controlled with control_bit_size qubits
    The input should be the quantum number register (number b), control qubits, classical registers (numbers a and N) 
    Control quibits are first control_bit_size 
    Quantum number b is the next bit_size quibits
    Classical number a is the first bit_size bits of classical register
    Classical number N is the next bit_size bits of classical register
    

    Parameter:
    __________
    bit_size, int, size of the binary numbers we add with this adder
    control_bit_size, int, size of the control qubit register
    
    Returns:
    __________
    qiskit.QuantumCircuit,instruction, a new instruction with adder modulo N
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    
    
    q_reg = qiskit.QuantumRegister(bit_size + ctrl_bit_size + 1) #this quantum register stores the added quantum number a
    cl_reg = qiskit.ClassicalRegister(2 * bit_size) #this classical register stores classical number b
    adder_mod_N = qiskit.QuantumCircuit(q_reg, cl_reg)

    inst = draper_adder_cl_instr(bit_size, ctrl_bit_size)
    qubits = []
    for i in range(0, ctrl_bit_size):
        qubits.append(q_reg[i])
    for i in range(bit_size):
        qubits.append(q_reg[i + ctrl_bit_size])
    cl_reg_a = []
    for i in range(0, bit_size):
        cl_reg_a.append(cl_reg[i])
    adder_mod_N.append(inst, qubits , cl_reg_a)

    cl_reg_N = []
    qubits = []
    for i in range(0, bit_size):
        qubits.append(q_reg[i + ctrl_bit_size])
    for i in range(0, bit_size):
        cl_reg_N.append(cl_reg[i + bit_size])
    inst = draper_subtraction_cl_instr(bit_size, 0)
    adder_mod_N.append(inst, qubits , cl_reg_N)

    q_anc = q_reg[bit_size + ctrl_bit_size]
    adder_mod_N.cx(qubits[0], q_anc)

    inst = draper_adder_cl_instr(bit_size, 1)
    qubits = [q_anc]
    for i in range(bit_size):
        qubits.append(q_reg[i + ctrl_bit_size])
    adder_mod_N.append(inst, qubits , cl_reg_N)

    inst = draper_subtraction_cl_instr(bit_size, ctrl_bit_size)
    qubits = []
    for i in range(0, ctrl_bit_size):
        qubits.append(q_reg[i])
    for i in range(bit_size):
        qubits.append(q_reg[i + ctrl_bit_size])   
    adder_mod_N.append(inst, qubits , cl_reg_a)

    adder_mod_N.x(q_reg[ctrl_bit_size])
    adder_mod_N.cx(q_reg[ctrl_bit_size], q_anc)
    adder_mod_N.x(q_reg[ctrl_bit_size])

    inst = draper_adder_cl_instr(bit_size, ctrl_bit_size)
    qubits = []
    for i in range(0, ctrl_bit_size):
        qubits.append(q_reg[i])
    for i in range(bit_size):
        qubits.append(q_reg[i + ctrl_bit_size])   
    adder_mod_N.append(inst, qubits , cl_reg_a)

    instr = adder_mod_N.to_instruction(None, 'add_mod_N \n ctrl: ' + str(ctrl_bit_size) + '\n bits: ' + str(bit_size))

    return instr

def ctrl_subtr_mod_N(bit_size, ctrl_bit_size):
    """the function implements the subtraction modulo N (N is classical number) which subtract a clasical number a from quantum number (b - a) mod N
    The algorithm follows the paper arXiv:quant-ph/0205095
    The function creates an instruction (single block) having bit_size + control_bit_size qubits and 2 * bit_size classical bits as input and as output
    Each number has bit_size bits
    The adder is controlled with control_bit_size qubits
    The input should be the quantum number register (number b), control qubits, classical registers (numbers a and N) 
    Control quibits are first control_bit_size 
    Quantum number b is the next bit_size quibits
    Classical number a is the first bit_size bits of classical register
    Classical number N is the next bit_size bits of classical register
    

    Parameter:
    __________
    bit_size, int, size of the binary numbers we add with this adder
    control_bit_size, int, size of the control qubit register
    
    Returns:
    __________
    qiskit.QuantumCircuit,instruction, a new instruction with adder modulo N
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    
    
    q_reg = qiskit.QuantumRegister(bit_size + ctrl_bit_size + 1) #this quantum register stores the added quantum number a
    cl_reg = qiskit.ClassicalRegister(2 * bit_size) #this classical register stores classical number b
    adder_mod_N = qiskit.QuantumCircuit(q_reg, cl_reg)

    inst = draper_subtraction_cl_instr(bit_size, ctrl_bit_size)
    qubits = []
    for i in range(0, ctrl_bit_size):
        qubits.append(q_reg[i])
    for i in range(bit_size):
        qubits.append(q_reg[i + ctrl_bit_size])
    cl_reg_a = []
    for i in range(0, bit_size):
        cl_reg_a.append(cl_reg[i])
    adder_mod_N.append(inst, qubits , cl_reg_a)

    cl_reg_N = []
    qubits = []
    for i in range(0, bit_size):
        qubits.append(q_reg[i + ctrl_bit_size])
    for i in range(0, bit_size):
        cl_reg_N.append(cl_reg[i + bit_size])
    inst = draper_subtraction_cl_instr(bit_size, 0)
    adder_mod_N.append(inst, qubits , cl_reg_N)

    q_anc = q_reg[bit_size + ctrl_bit_size]
    adder_mod_N.cx(qubits[0], q_anc)

    inst = draper_adder_cl_instr(bit_size, 1)
    qubits = [q_anc]
    for i in range(bit_size):
        qubits.append(q_reg[i + ctrl_bit_size])
    adder_mod_N.append(inst, qubits , cl_reg_N)

    inst = draper_adder_cl_instr(bit_size, ctrl_bit_size)
    qubits = []
    for i in range(0, ctrl_bit_size):
        qubits.append(q_reg[i])
    for i in range(bit_size):
        qubits.append(q_reg[i + ctrl_bit_size])   
    adder_mod_N.append(inst, qubits , cl_reg_a)

    adder_mod_N.x(q_reg[ctrl_bit_size])
    adder_mod_N.cx(q_reg[ctrl_bit_size], q_anc)
    adder_mod_N.x(q_reg[ctrl_bit_size])

    inst = draper_subtraction_cl_instr(bit_size, ctrl_bit_size)
    qubits = []
    for i in range(0, ctrl_bit_size):
        qubits.append(q_reg[i])
    for i in range(bit_size):
        qubits.append(q_reg[i + ctrl_bit_size])   
    adder_mod_N.append(inst, qubits , cl_reg_a)

    instr = adder_mod_N.to_instruction(None, 'subtr_mod_N \n ctrl: ' + str(ctrl_bit_size) + '\n bits: ' + str(bit_size))

    return instr

def ctrl_mult_mod_N(bit_size, cl_num_a, cl_num_N):
    """The function does (b+a*x)mod N operation, where b and x are quantum bit_size-bits numbers, a and N are classical numbers with n-bits size 
    The algorithm is taken from arXiv:quant-ph/0205095 
    The function creates an instruction (single block) having 4 * (bit_size+1) + 2 qubits and 2 * (bit_size+1) classical bits as input and as output
    Each number (a, b, x and N) has bit_size+1 bits. So, the size of the numbers used in the calculations is 1 bit bigger than the size of numbers used by the user of this function
    The multiplier is controlled with one qubit
    The input should be in the order from first to last:
    the control qubit, 
    the quantum numbers x,
    the quantum number b,
    the ancilla qubits for classical number a
    the ancilla qubits for classical number N
    additional ancilla qubit
    the classical register for number a
    the classical register for number b

    The function accepts clasical numbers a and N as parameters. Nevertheless they are converted into quantun bits inside the algorithm and then converted into a classical register.
    Quantum numbers x and b are should come from previous stages of the bigger algorithm.
    

    Parameter:
    __________
    bit_size, int, size of the binary numbers we add with this adder
    cl_num_a, list of 0 and 1, classical number a, with the first element of the list corresponding to the highest bit
    cl_num_N, list of 0 and 1, classical number N, with the first element of the list corresponding to the highest bitN
    
    Returns:
    __________
    qiskit.QuantumCircuit,instruction, a new instruction with the multiplier
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    bit_size = bit_size + 1
    cl_num_a1 = [0] + cl_num_a
    cl_num_N1 = [0] + cl_num_N
    
    q_reg = qiskit.QuantumRegister(4 * bit_size + 2) #this quantum register stores the numbers b and x
    cl_reg = qiskit.ClassicalRegister(2 * bit_size) #this classical register stores classical numbers a and N
    mult_mod_N = qiskit.QuantumCircuit(q_reg, cl_reg)

    #preparing classical number a
    qubits = []
    for i in range(bit_size):
        qubits.append(q_reg[2 * bit_size +1 + i])
    clbits = []
    for i in range(bit_size):
        clbits.append(cl_reg[i])
   
    mult_mod_N = qubit_binary_prepare(qubits, cl_num_a1, mult_mod_N)
    mult_mod_N = qubits_meas(qubits, clbits, mult_mod_N)
    mult_mod_N = qubit_binary_prepare(qubits, cl_num_a1, mult_mod_N)

    
    #preparing classical number N 
    qubits = []
    for i in range(bit_size):
        qubits.append(q_reg[3 * bit_size + 1 + i])
    clbits = []
    for i in range(bit_size):
        clbits.append(cl_reg[bit_size + i])
    mult_mod_N = qubit_binary_prepare(qubits, cl_num_N1, mult_mod_N)
    mult_mod_N = qubits_meas(qubits, clbits, mult_mod_N)
    mult_mod_N = qubit_binary_prepare(qubits, cl_num_N1, mult_mod_N)
    mult_mod_N.barrier()
    
    
    for i in range(bit_size-1):
        inst = ctrl_add_mod_N(bit_size, 2)
        qubits = [q_reg[0]]
        qubits.append(q_reg[bit_size - i])
        for j in range(bit_size):
            qubits.append(q_reg[bit_size + 1 + j])
        qubits.append(q_reg[4 * bit_size + 1])
        mult_mod_N.append(inst, qubits , cl_reg)
        if i < bit_size-2:
            qubits = []
            for j in range(bit_size):
                qubits.append(q_reg[2 * bit_size + 1 + j])
            #preparing classical number a * 2^i
            for j in range(bit_size - 1):
                cl_num_a1[j] = cl_num_a1[j + 1]
            cl_num_a1[bit_size - 1] = 0
            mult_mod_N = qubit_binary_prepare(qubits, cl_num_a1, mult_mod_N)
            clbits = []
            for i in range(bit_size):
                clbits.append(cl_reg[i])
            mult_mod_N = qubits_meas(qubits, clbits, mult_mod_N)
            mult_mod_N = qubit_binary_prepare(qubits, cl_num_a1, mult_mod_N)
            mult_mod_N.barrier()
      

    instr = mult_mod_N.to_instruction(None, 'amult_mod_N \n bits: ' + str(bit_size-1))

    return instr


def ctrl_mult_mod_N_s(bit_size, cl_num_a, cl_num_N):
    """The function does (b-a*x)mod N operation, where b and x are quantum bit_size-bits numbers, a and N are classical numbers with n-bits size 
    The algorithm is taken from arXiv:quant-ph/0205095 
    The function creates an instruction (single block) having 4 * (bit_size+1) + 2 qubits and 2 * (bit_size+1) classical bits as input and as output
    Each number (a, b, x and N) has bit_size+1 bits. So, the size of the numbers used in the calculations is 1 bit bigger than the size of numbers used by the user of this function
    The multiplier is controlled with one qubit
    The input should be in the order from first to last:
    the control qubit, 
    the quantum numbers x,
    the quentum number b,
    the ancilla qubits for classical number a
    the ancilla qubits for classical number N
    additional ancilla qubit
    the classical register for number a
    the classical register for number b

    The function accepts clasical numbers a and N as parameters. Nevertheless they are converted into quantun bits inside the algorithm and then converted into a classical register.
    Quantum numbers x and b  should come from previous stages of the bigger algorithm.
    

    Parameter:
    __________
    bit_size, int, size of the binary numbers we add with this adder
    cl_num_a, list of 0 and 1, classical number a, with the first element of the list corresponding to the highest bit
    cl_num_N, list of 0 and 1, classical number N, with the first element of the list corresponding to the highest bitN
    
    Returns:
    __________
    qiskit.QuantumCircuit,instruction, a new instruction with the multiplier
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """
    bit_size = bit_size + 1
    cl_num_a1 = [0] + cl_num_a
    cl_num_N1 = [0] + cl_num_N
    
    q_reg = qiskit.QuantumRegister(4 * bit_size + 2) #this quantum register stores the numbers b and x
    cl_reg = qiskit.ClassicalRegister(2 * bit_size) #this classical register stores classical numbers a and N
    mult_mod_N = qiskit.QuantumCircuit(q_reg, cl_reg)
    #preparing classical number a
    qubits = []
    for i in range(bit_size):
        qubits.append(q_reg[2 * bit_size + 1 + i])
    clbits = []
    for i in range(bit_size):
        clbits.append(cl_reg[i])
    mult_mod_N = qubit_binary_prepare(qubits, cl_num_a1, mult_mod_N)
    mult_mod_N = qubits_meas(qubits, clbits, mult_mod_N)
    mult_mod_N = qubit_binary_prepare(qubits, cl_num_a1, mult_mod_N)
    
    #preparing classical number N 
    qubits = []
    for i in range(bit_size):
        qubits.append(q_reg[3 * bit_size + 1 + i])
    clbits = []
    for i in range(bit_size):
        clbits.append(cl_reg[bit_size + i])
    mult_mod_N = qubit_binary_prepare(qubits, cl_num_N1, mult_mod_N)
    mult_mod_N = qubits_meas(qubits, clbits, mult_mod_N)
    mult_mod_N = qubit_binary_prepare(qubits, cl_num_N1, mult_mod_N)
    mult_mod_N.barrier()
    
    
    for i in range(bit_size-1):
        inst = ctrl_subtr_mod_N(bit_size, 2)
        qubits = [q_reg[0]]
        qubits.append(q_reg[bit_size - i])
        for j in range(bit_size):
            qubits.append(q_reg[bit_size + 1 + j])
        qubits.append(q_reg[4 * bit_size + 1])
        mult_mod_N.append(inst, qubits , cl_reg)
        if i < bit_size-2:
            qubits = []
            for j in range(bit_size):
                qubits.append(q_reg[2 * bit_size + 1 + j])
            #preparing classical number a * 2^i
            for j in range(bit_size - 1):
                cl_num_a1[j] = cl_num_a1[j + 1]
            cl_num_a1[bit_size - 1] = 0
            mult_mod_N = qubit_binary_prepare(qubits, cl_num_a1, mult_mod_N)
            clbits = []
            for i in range(bit_size):
                clbits.append(cl_reg[i])
            mult_mod_N = qubits_meas(qubits, clbits, mult_mod_N)
            mult_mod_N = qubit_binary_prepare(qubits, cl_num_a1, mult_mod_N)
            mult_mod_N.barrier()



    instr = mult_mod_N.to_instruction(None, 'amult_mod_N_s \n bits: ' + str(bit_size-1))

    return instr


def CQA_gate(bit_size, cl_num_a, cl_num_N):
    """The function does (a*x)mod N operation, where x is quantum bit_size-bits numbers, a and N are classical numbers with n-bits size 
    The algorithm is taken from arXiv:quant-ph/0205095 
    The difference comparing to function ctrl_mult_mod_N is that the output qubits are the same qubits used for input of x nuber. 
    All the ancilla qubits stay in state |0> after this gate.
    The function creates an instruction (single block) having 4 * (bit_size+1) + 2 qubits and 2 * (bit_size+1) classical bits as input and as output
    Each number (a, x and N) has bit_size+1 bits. 
    So, the size of the numbers used in the calculations is 1 bit bigger than the size of numbers used by the user of this function
    The multiplier is controlled with one qubit
    The input should be in the order from first to last:
    the control qubit, 
    the quantum numbers x,
    ancilla qubits (bit_size + 1 ),
    the ancilla qubits for classical number a
    the ancilla qubits for classical number N
    additional ancilla qubit
    the classical register for number a
    the classical register for number N

    The function accepts clasical numbers a and N as parameters. 
    Nevertheless they are converted into quantun bits inside the algorithm and then converted into a classical register.
    Quantum numbers x should come from previous stages of the bigger algorithm.
    

    Parameter:
    __________
    bit_size, int, size of the binary numbers we add with this adder
    cl_num_a, list of 0 and 1, classical number a, with the first element of the list corresponding to the highest bit
    cl_num_N, list of 0 and 1, classical number N, with the first element of the list corresponding to the highest bitN
    
    Returns:
    __________
    qiskit.QuantumCircuit,instruction, a new instruction with the multiplier
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    References:
    ___________
    My primitive study of the QFT and IQFT can be found at https://github.com/OlegUdalov/QC-qiskit-codes
    """

    bit_size = bit_size + 1
    cl_num_a = cl_num_a
    cl_num_N = cl_num_N
    
    q_reg = qiskit.QuantumRegister(4 * bit_size + 2) #this quantum register stores the numbers b and x
    cl_reg = qiskit.ClassicalRegister(2 * bit_size) #this classical register stores classical numbers a and N
    CQA_gate = qiskit.QuantumCircuit(q_reg, cl_reg)

    instr = ctrl_mult_mod_N(bit_size-1, cl_num_a, cl_num_N)
    CQA_gate.append(instr, q_reg, cl_reg)
    
    for i in range(bit_size):
        CQA_gate.swap(q_reg[i+1],q_reg[i+1+bit_size])

    a_int = cr.bin_str_2int(cl_num_a)
    N_int = cr.bin_str_2int(cl_num_N)
    a_inv_int = cr.euclids_inverse_mod(a_int, N_int)
    a_inv = cr.int_2_bin_str(a_inv_int, bit_size - 1)
    instr = ctrl_mult_mod_N_s(bit_size-1, a_inv, cl_num_N)
    CQA_gate.append(instr, q_reg, cl_reg)

    instr = CQA_gate.to_instruction(None, 'CQA_gate \n bits: ' + str(bit_size-1))
    
    return instr
    

def qubit_binary_prepare(q_reg, cl_reg, circ):
    """the function convert classical binary number into a 'quantum number'. If we have classical number 101, the function convert initial state 
    |0>|0>|0> into |1>|0>|1>. Application of the function to a superposition state will give wrong results. It is always applied to the initial state.
    
    If number of classical bits is smaller than number of qubits, the classical number is extended on the highest bits side.
    In the opposite case the classical number is cut on the highest bit side.

    Parameter:
    __________
    q_reg: qiskit.QuantumRegister, list of qubits that will be used modified (the first bit in the list represents the highest digit of quantum number)
    cl_reg: list, clsiical binary number (the first bit in the list represents the highest digit of classical number)
    circ: qiskit.QuantumCircuit, a quantum circuit to which the modification is made
    
    
    Returns:
    __________
    qiskit.QuantumCircuit, a new quantum circuit with modified qubits
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    """
    bit_size = len(q_reg) #size (number of digits) of the binary number
    if bit_size > len(cl_reg):
        diff = bit_size - len(cl_reg)
        ext = [0] * diff
        cl_reg = ext + cl_reg
        print('qubit_binary_prepare: Size of quantum number is bigger than that of classical, so the classical number extended by adding few highest bits')
        
    if bit_size < len(cl_reg):
        diff =  len(cl_reg) - bit_size
        del cl_reg[0 : diff]
        print('qubit_binary_prepare: Size of quantum number is smaller than that of classical, so the classical number is cut (highest bits)')

    for i in range(bit_size):
        if cl_reg[i] == 1:
            circ.x(q_reg[i])

    return circ


def qubits_meas(q_reg, cl_reg, circ):
    """the function adds measurement of qubits to the circuit q_reg to the classical bits cl_reg
    Number of qubits should be the same as number of classical registers. Othewise the circuit do nothing and returns a message
 
    
    Parameter:
    __________
    q_reg: qiskit.QuantumRegister, list of qubits that will be used measured
    cl_reg: qiskit.ClassicalRegister, clzssical registers where the measurement results will be stored
    circ: qiskit.QuantumCircuit, a quantum circuit to which the modification is made
    
    
    Returns:
    __________
    qiskit.QuantumCircuit, a new quantum circuit added measurements
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    
    """

    bit_size = len(q_reg) #size (number of digits) of the binary number
    cl_reg_size = len(cl_reg)
    if bit_size != cl_reg_size:
        print('qubits_meas: Size of quantum register is not the same as size of classical one, so the function does nothing')
        return circ
     
    for i in range(bit_size):
        circ.measure(q_reg[i], cl_reg[i])

    return circ

def plot_hyst_func(circ, shot_num):
    ''' the function performs quantum algorithm multiple times and plot a hystogram of measured state
    Note that the function can be used if only there are measuremenst in the circuit
    
    Parameter:
    __________
    circ: qiskit.QuantumCircuit, quantum circuit to run
    shot_num: int, number of shots of the algorithm
    
    
    Returns:
    __________
    plots a hystogram
        
    Notes:
    __________
    Qiskit 1.0 was used to create the function
    '''
    simulator = Aer.get_backend('qasm_simulator') 
    simulator = AerSimulator()
    circ = transpile(circ, simulator)
    result = simulator.run(circ,shots = shot_num).result()
    counts = result.get_counts(circ)
    return plot_histogram(counts, title = 'results')
    

