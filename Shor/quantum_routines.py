'''this file contains some quantum routines that I will use in Shor algorithm'''


import qiskit
import numpy as np

def QFTn(n, Psi_init, Measure):
    '''this function implements the QFT circuit'''
    QubitsNumber=n
    qubits1 = qiskit.QuantumRegister(QubitsNumber)
    if Measure=='yes':
        Readout=qiskit.ClassicalRegister(QubitsNumber)
        QFT = qiskit.QuantumCircuit(qubits1,Readout)
    else:
        QFT = qiskit.QuantumCircuit(qubits1)
    
    norm=np.linalg.norm(Psi_init)
    Psi_init=Psi_init/norm
    QFT.initialize(Psi_init)

    for i in range(0,QubitsNumber):
        QFT.h(qubits1[QubitsNumber-1-i])
        if i<QubitsNumber-1:
            for j in range(0,QubitsNumber-i-1):
                Rangle=PI/(math.pow(2,j+1))
                QFT.cp(Rangle,qubits1[QubitsNumber-1-i],qubits1[QubitsNumber-2-j-i])
    for i in range(0,QubitsNumber//2):
        QFT.swap(qubits1[i],qubits1[QubitsNumber-1-i])
    
    if Measure=='yes':
        QFT.measure(qubits1, Readout)
        return [QFT,qubits1,Readout]
    else:
        return [QFT,qubits1]