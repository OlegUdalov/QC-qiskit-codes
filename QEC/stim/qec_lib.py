# QEC functions for stim package 

import stim


def synd_x(qubits):
    '''
    The function creates X syndrome quantum circuit
    The circuit consists of a series of CX gates controlled by the same qubit - last qubit of the input qubit array
    All qubits beside the last one are the target qubits

    Input parametes:
    qubits, int array - set of qubit numbers

    Output:
    quantum circuit that can be added to the full QEC code
    '''
    qc = stim.Circuit()
    n_qubits = max(qubits)
    size = len(qubits)
    q = range(n_qubits + 1)
    if size > 2:
        for i_q in range(size - 1):
            qc.append_operation('CX',[q[qubits[size - 1]], q[qubits[i_q]]])
    else:
        qc.append_operation('CX',[q[qubits[size - 1]], q[qubits[0]]])
    return qc

def synd_z(qubits):
    '''
    The function creates Z syndrome quantum circuit
    The circuit consists of a series of CZ gates controlled by the same qubit - last qubit of the input qubit array
    All qubits beside the last one are the target qubits

    Input parametes:
    qubits, int array - set of qubit numbers

    Output:
    quantum circuit that can be added to the full QEC code
    '''
    qc = stim.Circuit()
    n_qubits = max(qubits)
    size = len(qubits)
    q = range(n_qubits + 1)
    if size > 2:
        for i_q in range(size - 1):
            qc.append_operation('CZ',[q[qubits[size - 1]], q[qubits[i_q]]])
    else:
        qc.append_operation('CZ',[q[qubits[size - 1]], q[qubits[0]]])
    return qc