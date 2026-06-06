OPENQASM 2.0;
include "qelib1.inc";
// steane_713_transversal_H: [[7,1,3]] Steane -- transversal logical H = H on all 7 qubits (self-dual CSS); realizes logical H
qreg q[7];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
