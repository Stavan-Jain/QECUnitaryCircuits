OPENQASM 2.0;
include "qelib1.inc";
// ghz5_fanout: 5-qubit GHZ state -- GHZ via H + CNOT fan-out from q0 (alternate construction to the ladder)
qreg q[5];
h q[0];
cx q[0],q[1];
cx q[0],q[2];
cx q[0],q[3];
cx q[0],q[4];
