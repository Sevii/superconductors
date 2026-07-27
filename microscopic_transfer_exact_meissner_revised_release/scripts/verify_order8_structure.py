#!/usr/bin/env python3
"""Certificate for the complete degree-eight source formula and obstruction.

The script verifies, on every product-AGP filling of the three-cell ring, that
for the single-flavor screened model

 R8 P_Z = b^2 U P_Z + b^4 V P_Z,
 U P_Z = -16 s^4/(dm^3 ds^2) (sum_e B_e^2) P_Z,
 V P_Z = dm^-3 M P_Z,
 M = sum_e B_e^4 + 1/2 sum_{e<f} [B_e,B_f]^* [B_e,B_f].

It also verifies the Kramers closed-shell replacement
 M -> M_Theta = 1/2 sum_e B_e^4
                  + 1/2 sum_{e<f} [B_e,B_f]^* [B_e,B_f],
the separate vanishing of the commutator-square matrix element for the
specified two-pair witness, and the exact amplitudes 8 sqrt(3) and 4 sqrt(3).
"""
from __future__ import annotations
import importlib.util
import math
from pathlib import Path
import sys

import numpy as np
import scipy
import scipy.linalg as la

ROOT=Path(__file__).parent
spec=importlib.util.spec_from_file_location('o8',ROOT/'audit_order8_all_fillings.py')
o8=importlib.util.module_from_spec(spec);sys.modules[spec.name]=o8;assert spec.loader;o8.__spec__=spec;spec.loader.exec_module(o8)


def opnorm(x): return 0.0 if x.size==0 else float(la.norm(x,2))


def raw(active,ctrl,dm,b,s,a):
    C2,A,Bq,Cq,r0=o8.build_operators(active,ctrl,b,s,a,dm)
    Z=active.Z
    A4=-o8.word_action(C2,(A,A,A,A),r0,Z)
    B2=-o8.word_action(C2,(Bq,Bq),r0,Z)
    C=+o8.word_action(C2,(Cq,),r0,Z)
    AAB=sum((o8.word_action(C2,w,r0,Z) for w in ((A,A,Bq),(A,Bq,A),(Bq,A,A))),start=np.zeros_like(Z))
    return A4,B2+C,AAB


def main():
    L=3;dm=10.0;ds=8.0;da=8.0;b=0.35;s=0.9;a=0.45
    ctrl=o8.control_space(L,dm,ds,da)
    errors={"router_source":0.0,"bridge_metric":0.0,"mixed":0.0,"completion":0.0}
    print('COMPLETE DEGREE-EIGHT STRUCTURE CERTIFICATE')
    print('='*78)
    print('numpy',np.__version__,'scipy',scipy.__version__)
    print(f'L={L}, dm={dm}, ds={ds}, da={da}, b={b}, s={s}, a={a}')
    print()

    witness=None
    for n in range(2*L+1):
        active=o8.active_sector(L,2*n)
        Z=active.Z
        A4,BC,AAB=raw(active,ctrl,dm,b,s,a)
        Bsq=sum((B@B for B in active.B),start=0*active.eye)
        BsqZ=Bsq@Z
        MZ=np.zeros_like(Z)
        localZ=np.zeros_like(Z)
        commZ=np.zeros_like(Z)
        for B in active.B:
            B4Z=(B@B@B@B)@Z
            localZ+=B4Z
            MZ+=B4Z
        for e in range(L):
            for f in range(e+1,L):
                comm=active.B[e]@active.B[f]-active.B[f]@active.B[e]
                term=0.5*(comm.getH()@comm)@Z
                commZ+=term
                MZ+=term
        expected_A4=-(16*b*b*s**4/(dm**3*ds**2))*BsqZ
        expected_BC=(b**4/dm**3)*MZ
        completion= A4+BC + (16*b*b*s**4/(dm**3*ds**2))*BsqZ - (b**4/dm**3)*MZ
        ea=opnorm(A4-expected_A4); ev=opnorm(BC-expected_BC); em=opnorm(AAB); ec=opnorm(completion)
        errors['router_source']=max(errors['router_source'],ea)
        errors['bridge_metric']=max(errors['bridge_metric'],ev)
        errors['mixed']=max(errors['mixed'],em)
        errors['completion']=max(errors['completion'],ec)
        print(f'n={n}: err(A4 source)={ea:.3e}, err(bridge metric)={ev:.3e}, mixed={em:.3e}, completed source={ec:.3e}')

        if n==2:
            # Input |0,2>; output d^+_{1,0}d^+_{1,1}|0> in block 1 (zero-based block 0).
            col=active.compositions.index((0,2))
            state=0
            for x in (0,1):
                state |= 1<<o8.mode(L,0,x,0)
                state |= 1<<o8.mode(L,0,x,1)
            row=active.basis.index(state)
            MthetaZ=0.5*localZ+commZ
            witness={
                'B2':BsqZ[row,col],
                'local':localZ[row,col],
                'comm':commZ[row,col],
                'M':MZ[row,col],
                'Mtheta':MthetaZ[row,col],
            }

    print()
    assert witness is not None
    print('Two-pair all-filling obstruction witness')
    print(f"  <Phi|sum B_e^2|0,2>       = {witness['B2'].real:+.12e}")
    print(f"  <Phi|sum B_e^4|0,2>       = {witness['local'].real:+.12e}  (expected 8 sqrt(3))")
    print(f"  <Phi|commutator sum|0,2>  = {witness['comm'].real:+.12e}  (expected 0)")
    print(f"  <Phi|M|0,2>               = {witness['M'].real:+.12e}  (expected 8 sqrt(3))")
    print(f"  <Phi|M_Theta|0,2>         = {witness['Mtheta'].real:+.12e}  (expected 4 sqrt(3))")
    ew=max(abs(witness['B2']),abs(witness['local']-8*math.sqrt(3)),abs(witness['comm']),abs(witness['M']-8*math.sqrt(3)),abs(witness['Mtheta']-4*math.sqrt(3)))
    print(f'  witness error             = {ew:.3e}')
    print()
    for k,v in errors.items(): print(f'max {k:16s} error = {v:.12e}')
    checks=[max(errors.values())<2e-9,ew<2e-9,abs(witness['M'])>1]
    print('Conclusion: neither the single-flavor nor closed-shell Kramers unmodified family can preserve the complete product-AGP zero manifold at weighted degree eight for b != 0.')
    print('OVERALL:','PASS' if all(checks) else 'FAIL')
    return 0 if all(checks) else 1

if __name__=='__main__': raise SystemExit(main())
