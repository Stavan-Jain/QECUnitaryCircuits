OPENQASM 2.0;
include "qelib1.inc";
// rep3_phaseflip_encode: [[3,1,1]] phase-flip repetition -- encoder: CNOT fan-out conjugated by Hadamards
qreg q[3];
cx q[0],q[1];
cx q[0],q[2];
h q[0];
h q[1];
h q[2];
