OPENQASM 2.0;
include "qelib1.inc";
// steane_713_transversal_CNOT: [[7,1,3]] Steane x2 -- transversal logical CNOT between two Steane blocks (qubits 0-6 control, 7-13 target)
qreg q[14];
cx q[0],q[7];
cx q[1],q[8];
cx q[2],q[9];
cx q[3],q[10];
cx q[4],q[11];
cx q[5],q[12];
cx q[6],q[13];
