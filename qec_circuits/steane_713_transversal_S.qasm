OPENQASM 2.0;
include "qelib1.inc";
// steane_713_transversal_S: [[7,1,3]] Steane -- transversal S(dg) layer; realizes logical S (phase gate) on the Steane code
qreg q[7];
sdg q[0];
sdg q[1];
sdg q[2];
sdg q[3];
sdg q[4];
sdg q[5];
sdg q[6];
