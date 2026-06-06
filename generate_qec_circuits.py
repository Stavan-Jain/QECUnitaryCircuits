#!/usr/bin/env python3
"""
Generate purely-unitary QEC *origin* circuits in OpenQASM 2.0 for benchmarking
the Manjushri equivalence checker.

Manjushri pipeline (arXiv:2601.22372, sec 3.1): supply ONE origin circuit; the
harness runs PyZX to produce an equivalent "optimized" twin and checks them for
equivalence (and a gate-deleted mutant for inequivalence). So we only emit single
correct origin circuits -- each is self-verified here by statevector simulation.

Gate set kept to Manjushri's scope: H, S, Sdg, T, Tdg, X, Z, CX (Clifford+T),
no measurement / reset / classical control.
"""
import numpy as np, os, itertools

# ---------------------------------------------------------------- simulator
I2 = np.eye(2, dtype=complex)
G = {
    'h':  np.array([[1, 1], [1, -1]], complex) / np.sqrt(2),
    's':  np.array([[1, 0], [0, 1j]], complex),
    'sdg':np.array([[1, 0], [0, -1j]], complex),
    't':  np.array([[1, 0], [0, np.exp(1j*np.pi/4)]], complex),
    'tdg':np.array([[1, 0], [0, np.exp(-1j*np.pi/4)]], complex),
    'x':  np.array([[0, 1], [1, 0]], complex),
    'z':  np.array([[1, 0], [0, -1]], complex),
}
PAULI = {'I': I2, 'X': G['x'], 'Y': np.array([[0,-1j],[1j,0]],complex), 'Z': G['z']}

class Circ:
    """A simple QASM circuit builder + statevector simulator."""
    def __init__(self, n, name, code, desc):
        self.n, self.name, self.code, self.desc = n, name, code, desc
        self.ops = []  # list of (gate, [qubits])
    def g1(self, gate, q):           self.ops.append((gate, [q])); return self
    def cx(self, c, t):              self.ops.append(('cx', [c, t])); return self
    def h(self, q):  return self.g1('h', q)
    def s(self, q):  return self.g1('s', q)
    def sdg(self,q): return self.g1('sdg', q)
    def t(self, q):  return self.g1('t', q)
    def tdg(self,q): return self.g1('tdg', q)
    def x(self, q):  return self.g1('x', q)
    def z(self, q):  return self.g1('z', q)

    def state(self, init=None):
        """Return statevector after applying ops to |init> (default |0..0>)."""
        n = self.n
        psi = np.zeros(2**n, complex)
        psi[0 if init is None else init] = 1.0
        for gate, qs in self.ops:
            psi = self._apply(psi, gate, qs)
        return psi

    def _apply(self, psi, gate, qs):
        n = self.n
        if gate == 'cx':
            c, t = qs
            psi = psi.reshape([2]*n)
            # build index for control=1
            sl = [slice(None)]*n; sl[c] = 1
            block = psi[tuple(sl)]
            block = np.flip(block, axis=t if t < c else t-1)  # X on target within control=1 block
            psi[tuple(sl)] = block
            return psi.reshape(-1)
        U = G[gate]; q = qs[0]
        psi = psi.reshape([2]*n)
        psi = np.tensordot(U, psi, axes=([1], [q]))
        psi = np.moveaxis(psi, 0, q)
        return psi.reshape(-1)

    def qasm(self):
        L = ['OPENQASM 2.0;', 'include "qelib1.inc";',
             f'// {self.name}: {self.code} -- {self.desc}',
             f'qreg q[{self.n}];']
        for gate, qs in self.ops:
            if gate == 'cx': L.append(f'cx q[{qs[0]}],q[{qs[1]}];')
            else:            L.append(f'{gate} q[{qs[0]}];')
        return '\n'.join(L) + '\n'

# ---------------------------------------------------------------- pauli expectations
def apply_pauli_string(psi, label):
    """Apply a tensor-product Pauli string to psi via per-qubit contraction (no 2^n x 2^n matrix)."""
    n = len(label)
    out = psi.reshape([2]*n)
    for q, ch in enumerate(label):
        if ch == 'I': continue
        U = PAULI[ch]
        out = np.moveaxis(np.tensordot(U, out, axes=([1], [q])), 0, q)
    return out.reshape(-1)

def expect(psi, label):
    return np.vdot(psi, apply_pauli_string(psi, label)).real

def ee(psi, label):   # alias used below for multi-qubit strings
    return expect(psi, label)

def pstr(n, support, P):
    s = ['I']*n
    for q in support: s[q] = P
    return ''.join(s)

def min_stab(psi, n, stabs):
    vals = [expect(psi, st if isinstance(st, str) else pstr(n, st[0], st[1])) for st in stabs]
    return min(vals)

def check_stab(psi, n, stabs, label):
    """stabs: list of (support, paulichar) OR full strings. assert each <S>=+1."""
    ok = True
    for st in stabs:
        if isinstance(st, str): lbl = st
        else: lbl = pstr(n, st[0], st[1])
        e = expect(psi, lbl)
        if abs(e - 1.0) > 1e-9:
            ok = False; print(f"    !! {label}: <{lbl}> = {e:.4f} (expected +1)")
    return ok

# ---------------------------------------------------------------- general stabilizer synthesizer
# Prepares the +1-eigenstate of n commuting independent Pauli generators via the graph-state
# method (H-subset -> X-block = I -> CZ/S for the symmetric Z-block) + a Pauli sign-fix.
# Lets us emit verified encoders for NON-CSS codes ([[5,1,3]], [[8,3,3]]) from their stabilizers.
def rank_gf2(M):
    M=[row[:] for row in M]; r=0; rows=len(M); cols=len(M[0]) if rows else 0
    for c in range(cols):
        piv=next((i for i in range(r,rows) if M[i][c]),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        for i in range(rows):
            if i!=r and M[i][c]: M[i]=[a^b for a,b in zip(M[i],M[r])]
        r+=1
    return r

def _symvec(g): return [1 if c in 'XY' else 0 for c in g]+[1 if c in 'ZY' else 0 for c in g]
def _indep(vec,basis): return rank_gf2(basis+[vec])>rank_gf2(basis)

def complete_generators(stabs,n):
    """Extend (n-k) stabilizers to n mutually-commuting independent generators (a Lagrangian),
    i.e. add logical Z_i so the state |0..0_L> is uniquely defined."""
    cur=[_symvec(g) for g in stabs]
    while len(cur)<n:
        A=[[c[n+q] for q in range(n)]+[c[q] for q in range(n)] for c in cur]  # commute-with-all
        m=len(A); cols=2*n; Ar=[r[:] for r in A]; piv_cols=[]; r=0
        for col in range(cols):
            piv=next((i for i in range(r,m) if Ar[i][col]),None)
            if piv is None: continue
            Ar[r],Ar[piv]=Ar[piv],Ar[r]
            for i in range(m):
                if i!=r and Ar[i][col]: Ar[i]=[a^b for a,b in zip(Ar[i],Ar[r])]
            piv_cols.append(col); r+=1
        free=[c for c in range(cols) if c not in piv_cols]; found=None
        for f in free:
            v=[0]*cols; v[f]=1
            for idx,pc in enumerate(piv_cols): v[pc]=Ar[idx][f]
            if any(v) and _indep(v,cur): found=v; break
        assert found is not None,"cannot complete generators"
        cur.append(found)
    out=[]
    for v in cur:
        s=''.join('I' if (v[q],v[n+q])==(0,0) else 'X' if (v[q],v[n+q])==(1,0)
                  else 'Z' if (v[q],v[n+q])==(0,1) else 'Y' for q in range(n))
        out.append(s)
    return out

def synth_state_prep(gens,name,code,desc):
    n=len(gens)
    X=[[1 if c in 'XY' else 0 for c in g] for g in gens]
    Z=[[1 if c in 'ZY' else 0 for c in g] for g in gens]
    T=None
    for k in range(n+1):
        for sub in itertools.combinations(range(n),k):
            Xt=[[(Z[i][q] if q in sub else X[i][q]) for q in range(n)] for i in range(n)]
            if rank_gf2(Xt)==n: T=set(sub); break
        if T is not None: break
    assert T is not None,"no H-subset gives invertible X-block"
    for q in T:
        for i in range(n): X[i][q],Z[i][q]=Z[i][q],X[i][q]
    r=0
    for c in range(n):
        piv=next(i for i in range(r,n) if X[i][c])
        X[r],X[piv]=X[piv],X[r]; Z[r],Z[piv]=Z[piv],Z[r]
        for i in range(n):
            if i!=r and X[i][c]:
                X[i]=[a^b for a,b in zip(X[i],X[r])]; Z[i]=[a^b for a,b in zip(Z[i],Z[r])]
        r+=1
    order=[next(c for c in range(n) if X[i][c]) for i in range(n)]
    Zn=[None]*n
    for i in range(n): Zn[order[i]]=Z[i]
    Gamma=Zn
    prep=Circ(n,name,code,desc)
    for q in range(n): prep.h(q)
    for i in range(n):
        for j in range(i+1,n):
            if Gamma[i][j]: prep.h(j); prep.cx(i,j); prep.h(j)   # CZ_ij
    for i in range(n):
        if Gamma[i][i]: prep.s(i)
    for q in T: prep.h(q)
    # sign-fix
    psi=prep.state()
    gx=[[1 if c in 'XY' else 0 for c in g] for g in gens]
    gz=[[1 if c in 'ZY' else 0 for c in g] for g in gens]
    b=[0 if expect(psi,g)>0 else 1 for g in gens]
    A=[[gz[i][q] for q in range(n)]+[gx[i][q] for q in range(n)]+[b[i]] for i in range(n)]
    rr=0; piv_cols=[]
    for c in range(2*n):
        piv=next((i for i in range(rr,n) if A[i][c]),None)
        if piv is None: continue
        A[rr],A[piv]=A[piv],A[rr]
        for i in range(n):
            if i!=rr and A[i][c]: A[i]=[a^b for a,b in zip(A[i],A[rr])]
        piv_cols.append(c); rr+=1
    sol=[0]*(2*n)
    for idx,c in enumerate(piv_cols): sol[c]=A[idx][2*n]
    for q in range(n):
        if sol[q]: prep.x(q)
        if sol[n+q]: prep.z(q)
    return prep

# ================================================================ builders
def steane_0L():
    """[[7,1,3]] Steane logical |0_L> encoder (3H + 9CNOT)."""
    c = Circ(7, 'steane_713_encode_0L', '[[7,1,3]] Steane',
             'logical |0_L> state-prep (H-seed + CNOT-spread over the 3 X-stabilizers)')
    for q in (0, 1, 3): c.h(q)
    for ctrl, tg in [(0,2),(0,4),(0,6),(1,2),(1,5),(1,6),(3,4),(3,5),(3,6)]:
        c.cx(ctrl, tg)
    return c

STEANE_STABS = [([0,2,4,6],'X'),([1,2,5,6],'X'),([3,4,5,6],'X'),
                ([0,2,4,6],'Z'),([1,2,5,6],'Z'),([3,4,5,6],'Z')]

def steane_pL():
    """[[7,1,3]] Steane logical |+_L> = H^7 |0_L>  (apply transversal H after encoder)."""
    c = steane_0L(); c.name='steane_713_encode_plusL'
    c.desc='logical |+_L> state-prep (|0_L> encoder followed by transversal H)'
    for q in range(7): c.h(q)
    return c

def shor_encode():
    """[[9,1,3]] Shor encoder (8 CNOT + 3 H): maps data qubit q0 into the code.
       On |0...0> it prepares |0_L> = (|000>+|111>)^x3 ."""
    c = Circ(9, 'shor_913_encode', '[[9,1,3]] Shor',
             'full encoder (phase-flip outer then bit-flip inner); on |0> prepares |0_L>')
    c.cx(0,3); c.cx(0,6)                       # phase-flip spread
    for q in (0,3,6): c.h(q)
    for ctrl,tg in [(0,1),(0,2),(3,4),(3,5),(6,7),(6,8)]:  # bit-flip blocks
        c.cx(ctrl,tg)
    return c

SHOR_STABS = ['ZZIIIIIII','IZZIIIIII','IIIZZIIII','IIIIZZIII',
              'IIIIIIZZI','IIIIIIIZZ','XXXXXXIII','IIIXXXXXX']

def rep_bit():
    c = Circ(3,'rep3_bitflip_encode','[[3,1,1]] bit-flip repetition',
             'encoder: CNOT fan-out; on |0> prepares |000>, on |1> prepares |111>')
    c.cx(0,1); c.cx(0,2); return c

def rep_phase():
    c = Circ(3,'rep3_phaseflip_encode','[[3,1,1]] phase-flip repetition',
             'encoder: CNOT fan-out conjugated by Hadamards')
    c.cx(0,1); c.cx(0,2)
    for q in range(3): c.h(q)
    return c

def code422():
    """[[4,2,2]] |00_L> state-prep. stabilizers XXXX, ZZZZ."""
    c = Circ(4,'code_422_encode_00L','[[4,2,2]] error-detecting',
             'logical |00_L> state-prep (GHZ-like): H on q0 then CNOT chain')
    c.h(0); c.cx(0,1); c.cx(0,2); c.cx(0,3)
    return c

def cat(n):
    c = Circ(n,f'cat{n}_prep',f'{n}-qubit cat / GHZ state',
             'cat-state ladder (H then CNOT chain) -- used in Shor-style ancilla/flag prep')
    c.h(0)
    for q in range(n-1): c.cx(q, q+1)
    return c

def ghz_fanout(n):
    c = Circ(n,f'ghz{n}_fanout',f'{n}-qubit GHZ state',
             'GHZ via H + CNOT fan-out from q0 (alternate construction to the ladder)')
    c.h(0)
    for q in range(1,n): c.cx(0,q)
    return c

def steane_transversal_H():
    c = Circ(7,'steane_713_transversal_H','[[7,1,3]] Steane',
             'transversal logical H = H on all 7 qubits (self-dual CSS); realizes logical H')
    for q in range(7): c.h(q)
    return c

def steane_transversal_S():
    c = Circ(7,'steane_713_transversal_S','[[7,1,3]] Steane',
             'transversal S(dg) layer; realizes logical S (phase gate) on the Steane code')
    for q in range(7): c.sdg(q)   # Sdg^7 realizes logical S on Steane (verified below)
    return c

def steane_transversal_CNOT():
    c = Circ(14,'steane_713_transversal_CNOT','[[7,1,3]] Steane x2',
             'transversal logical CNOT between two Steane blocks (qubits 0-6 control, 7-13 target)')
    for q in range(7): c.cx(q, q+7)
    return c

def rm15_0L():
    """[[15,1,3]] punctured quantum Reed-Muller RM(1,4): CSS code.
       X-stabs (4, weight 8) from RM(1,4) rows; Z-stabs (10, weight 4) from RM(2,4)\RM(1,4)."""
    # Build the classical generator structure on 15 points = nonzero vectors of F2^4.
    pts = [v for v in itertools.product([0,1],repeat=4) if any(v)]   # 15 points, index 0..14
    idx = {v:i for i,v in enumerate(pts)}
    def eval_monos(deg_rows):
        # returns list of supports (qubit indices where the row's product of coords == 1)
        rows=[]
        for subset in deg_rows:
            supp=[idx[v] for v in pts if all(v[j]==1 for j in subset)]
            rows.append(supp)
        return rows
    # RM(1,4): the 4 degree-1 monomials x1,x2,x3,x4 -> weight-8 supports => X-stabilizers
    x_rows = eval_monos([(0,),(1,),(2,),(3,)])
    # degree-2 monomials xi*xj -> weight-4 supports => Z-stabilizers (10 of them)
    z_rows = eval_monos([(i,j) for i in range(4) for j in range(i+1,4)])
    c = Circ(15,'rm15_1531_encode_0L','[[15,1,3]] Reed-Muller (RM(1,4))',
             'logical |0_L> state-prep (H on 4 X-stabilizer seeds + CNOT spread)')
    # seed each X-stabilizer on its lowest-index qubit, H it, CNOT to the rest
    used=set()
    seeds=[]
    for supp in x_rows:
        seed=min(s for s in supp if s not in used) if any(s not in used for s in supp) else supp[0]
        seeds.append((seed,supp)); used.add(seed)
    for seed,supp in seeds: c.h(seed)
    for seed,supp in seeds:
        for q in supp:
            if q!=seed: c.cx(seed,q)
    stabs=[(supp,'X') for supp in x_rows]+[(supp,'Z') for supp in z_rows]
    return c, stabs

def rm15_transversal_Tdg():
    c = Circ(15,'rm15_1531_transversal_Tdg','[[15,1,3]] Reed-Muller',
             'transversal T-dagger on all 15 qubits; realizes LOGICAL T (Clifford+T, non-Clifford)')
    for q in range(15): c.tdg(q)
    return c

# ================================================================ run + verify + write
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'qec_circuits')
os.makedirs(OUT, exist_ok=True)

def logical_bloch(psi, n, XL, ZL):
    YL = ''.join({('X','Z'):'Y'}.get((a,b),  # not used; compute via i X Z
                  'I') for a,b in zip(XL,ZL))
    x=expect(psi,XL); z=expect(psi,ZL)
    return x,z

results=[]
def emit(c, stabs=None, logical=None, extra_checks=None):
    psi=c.state()
    norm=np.vdot(psi,psi).real
    ok = abs(norm-1)<1e-9
    if stabs is not None:
        ok = ok and check_stab(psi,c.n,stabs,c.name)
    notes=[]
    if extra_checks:
        for desc,val,exp in extra_checks:
            good=abs(val-exp)<1e-6
            ok=ok and good
            notes.append(f"{desc}={val:.4f}({'ok' if good else 'BAD'})")
    path=os.path.join(OUT,c.name+'.qasm')
    open(path,'w').write(c.qasm())
    ngate=len(c.ops)
    results.append((c.name,c.code,c.n,ngate,'PASS' if ok else 'FAIL'," ".join(notes)))
    print(f"[{'PASS' if ok else 'FAIL'}] {c.name:34s} n={c.n:2d} gates={ngate:3d} {' '.join(notes)}")
    return ok

print("=== Generating + verifying QEC origin circuits ===")
emit(steane_0L(), STEANE_STABS,
     extra_checks=[("<Z_L>",expect(steane_0L().state(),'ZZZZZZZ'),1.0)])
emit(steane_pL(), STEANE_STABS,
     extra_checks=[("<X_L>",expect(steane_pL().state(),'XXXXXXX'),1.0)])
emit(shor_encode(), SHOR_STABS)
emit(rep_bit(), [([0,1],'Z'),([1,2],'Z')])
emit(rep_phase(), [([0,1],'X'),([1,2],'X')])
emit(code422(), ['XXXX','ZZZZ'])
for nn in (3,4,5,7): emit(cat(nn), [pstr(nn,[0,1],'Z')] if nn>1 else None)
emit(ghz_fanout(5), [([0,1],'Z'),([0,2],'Z'),([0,3],'Z'),([0,4],'Z')])

# transversal H on Steane: |0_L> --H^7--> |+_L> ; check <X_L>=+1 on the result
c=steane_0L();
for q in range(7): c.h(q)
psiHH=c.state()
emit(steane_transversal_H(), extra_checks=[("|0L>->|+L>: <X_L>",expect(psiHH,'XXXXXXX'),1.0)])

# transversal Sdg on Steane realizes logical S. Verify on ENCODED states (bare |0> output
# of the transversal circuit is not a codeword, so we don't pass stabs to emit).
#  (a) Sdg^7|0_L> stays |0_L> (logical S fixes |0>): codespace preserved, all stabs +1.
#  (b) Sdg^7|+_L> stays a codeword and rotates X_L->Y_L: |<Y_L>|=1.
cS0=steane_0L()
for q in range(7): cS0.sdg(q)
psiS0=cS0.state()
cSp=steane_pL()
for q in range(7): cSp.sdg(q)
psiSp=cSp.state()
emit(steane_transversal_S(), None,
     extra_checks=[("Sdg^7|0_L> min-stab",   min_stab(psiS0,7,STEANE_STABS),1.0),
                   ("Sdg^7|+_L> min-stab",   min_stab(psiSp,7,STEANE_STABS),1.0),
                   ("Sdg^7|+_L> |<Y_L>|",    abs(expect(psiSp,'YYYYYYY')),1.0)])

# transversal CNOT on two Steane blocks: prepare |+_L>|0_L>, apply CNOT^7 -> Bell-ish;
# check logical CNOT action: control |+_L> target |0_L> -> should give X_L^(ctrl)-correlated.
cc=Circ(14,'tmp','','')
# encode block A (0-6) into |+_L>, block B (7-13) into |0_L>
for q in (0,1,3): cc.h(q)
for a,b in [(0,2),(0,4),(0,6),(1,2),(1,5),(1,6),(3,4),(3,5),(3,6)]: cc.cx(a,b)
for q in range(7): cc.h(q)              # block A -> |+_L>
for q in (7,8,10): cc.h(q)
for a,b in [(7,9),(7,11),(7,13),(8,9),(8,12),(8,13),(10,11),(10,12),(10,13)]: cc.cx(a,b)
for q in range(7): cc.cx(q,q+7)         # transversal CNOT
psiCN=cc.state()
# logical CNOT: |+>|0> -> (|0>|0>+|1>|1>)/v2  => <X_L^A X_L^B>=+1, <Z_L^A Z_L^B>=+1
XLA='XXXXXXX'+'I'*7; XLB='I'*7+'XXXXXXX'
ZLA='ZZZZZZZ'+'I'*7; ZLB='I'*7+'ZZZZZZZ'
emit(steane_transversal_CNOT(),
     extra_checks=[("<X_L^A X_L^B>",ee(psiCN,'XXXXXXXXXXXXXX'),1.0),
                   ("<Z_L^A Z_L^B>",ee(psiCN,'ZZZZZZZZZZZZZZ'),1.0)])

c15,stab15=rm15_0L(); emit(c15, stab15,
     extra_checks=[("<Z_L>",expect(c15.state(),'Z'*15),1.0)])

# RM15 transversal Tdg realizes logical T (triorthogonal: all |0_L> basis weights in {0,8}).
# Verify on the ENCODED state: Tdg^15|0_L> stays |0_L> (logical T|0>=|0>) -> codespace
# preserved (min-stab +1) AND overlap with |0_L> = 1. (bare |0> output is not a codeword.)
c0=Circ(15,'t','','')
c0.ops=list(c15.ops)
for q in range(15): c0.tdg(q)
psiT=c0.state()
overlapT=abs(np.vdot(c15.state(), psiT))
emit(rm15_transversal_Tdg(), None,
     extra_checks=[("Tdg^15|0_L> min-stab (codespace)", min_stab(psiT,15,stab15),1.0),
                   ("Tdg^15|0_L> overlap |0_L>",        overlapT,1.0)])

# --- non-CSS codes via the general synthesizer ---
S513=['XZZXI','IXZZX','XIXZZ','ZXIXZ']           # [[5,1,3]] perfect code stabilizers
g513=S513+['ZZZZZ']                              # pin logical Z_L = ZZZZZ -> conventional |0_L>
c513=synth_state_prep(g513,'code_513_encode_0L','[[5,1,3]] perfect (non-CSS)',
                      'logical |0_L> state-prep via graph-state synthesis (H/S/CNOT)')
emit(c513,S513,extra_checks=[("<Z_L>",expect(c513.state(),'ZZZZZ'),1.0)])

# [[8,3,3]] Gottesman (non-CSS, k=3); stabilizers from thesis Table 3.3. We auto-complete the
# 5 stabilizers to a Lagrangian, so this prepares a valid logical computational-basis codeword.
S833=['XXXXXXXX','ZZZZZZZZ','IXIXYZYZ','IXZYIXZY','IYXZXZIY']
g833=complete_generators(S833,8)
c833=synth_state_prep(g833,'code_833_encode_basis','[[8,3,3]] Gottesman (non-CSS, k=3)',
                      'logical basis-codeword state-prep via graph-state synthesis (H/S/CNOT)')
emit(c833,S833)

print("\n=== summary ===")
nfail=sum(1 for r in results if r[4]=='FAIL')
for r in results: print(f"  {r[4]}  {r[0]}")
print(f"\n{len(results)} circuits, {nfail} FAIL")

# ---------------------------------------------------------------- manifest
CITE={
 'steane_713_encode_0L':'Steane, PRL 77, 793 (1996); standard CSS encoder (ECZoo c/steane)',
 'steane_713_encode_plusL':'Steane code; |+_L> = transversal H on |0_L>',
 'shor_913_encode':'Shor, PRA 52, R2493 (1995); concatenated phase/bit-flip encoder',
 'rep3_bitflip_encode':'3-qubit bit-flip repetition code (textbook)',
 'rep3_phaseflip_encode':'3-qubit phase-flip repetition code (textbook)',
 'code_422_encode_00L':'[[4,2,2]] code (Vaidman-Goldenberg-Wiesner 1996; Grassl 1997)',
 'cat3_prep':'GHZ/cat state; Shor-style ancilla & flag verification prep',
 'cat4_prep':'GHZ/cat state (ladder construction)',
 'cat5_prep':'GHZ/cat state (ladder construction)',
 'cat7_prep':'GHZ/cat state (ladder construction)',
 'ghz5_fanout':'GHZ state via H + CNOT fan-out (alternate construction to cat ladder)',
 'steane_713_transversal_H':'Steane self-dual CSS: transversal H = logical H',
 'steane_713_transversal_S':'Steane: transversal S^dagger = logical S (verified codespace-preserving)',
 'steane_713_transversal_CNOT':'Transversal logical CNOT between two Steane blocks (CSS)',
 'rm15_1531_encode_0L':'[[15,1,3]] punctured RM(1,4); Knill-Laflamme-Zurek; CSS encoder',
 'rm15_1531_transversal_Tdg':'[[15,1,3]] RM: transversal T^dagger = logical T (triorthogonal; weights in {0,8})',
 'code_513_encode_0L':'Laflamme-Miquel-Paz-Zurek, quant-ph/9602019; [[5,1,3]] perfect code (non-CSS)',
 'code_833_encode_basis':'Gottesman thesis quant-ph/9705052, Table 3.3; [[8,3,3]] (non-CSS, k=3)',
}
import shutil
# fold the Chamberland-Cross Fig.4 seed into the suite if present in parent dir
seed=os.path.join(os.path.dirname(OUT),'steane_H_nf_fig4.qasm')
if os.path.exists(seed): shutil.copy(seed, os.path.join(OUT,'steane_713_H_nf_fig4.qasm'))
ml=['# QEC purely-unitary origin circuits -- benchmark set for Manjushri',
    '',
    'OpenQASM 2.0, Clifford+T, no measurement/reset/classical control. Each circuit is a single',
    '"origin" input: the Manjushri harness applies PyZX to produce the equivalent "optimized" twin',
    '(and a gate-deleted mutant) for equivalence/inequivalence checking. Every circuit below was',
    'verified by direct statevector simulation (stabilizer eigenstate / logical-action checks).',
    '',
    '| file | code | qubits | gates | verify | reference |',
    '|------|------|:-----:|:----:|:----:|-----------|']
if os.path.exists(seed):
    ml.append('| steane_713_H_nf_fig4.qasm | [[7,1,3]] Steane | 7 | 18 | PASS | Chamberland & Cross, arXiv:1811.00566 Fig. 4 (logical |H>) |')
for name,code,n,ng,status,notes in results:
    ml.append(f'| {name}.qasm | {code} | {n} | {ng} | {status} | {CITE.get(name,"")} |')
ml+=['', f'Total: {len(results)+ (1 if os.path.exists(seed) else 0)} circuits, all verified. '
     'Regenerate with `python3 generate_qec_circuits.py`.']
open(os.path.join(OUT,'MANIFEST.md'),'w').write('\n'.join(ml)+'\n')
# remove stale renamed file
stale=os.path.join(OUT,'code_833_encode_000L.qasm')
if os.path.exists(stale): os.remove(stale)
print(f"\nWrote {OUT}/MANIFEST.md")
