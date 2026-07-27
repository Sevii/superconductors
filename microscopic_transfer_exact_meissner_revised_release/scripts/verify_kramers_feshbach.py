#!/usr/bin/env python3
"""Independent exact-Feshbach check of the Kramers bridge coefficient."""
from __future__ import annotations
import itertools
import math
import sys
import numpy as np
import scipy
import scipy.linalg as la

RNG=np.random.default_rng(20260727)


def herm(x): return (x+x.conj().T)/2

def opnorm(x): return float(la.norm(x,2))


def control_x(nchan,site):
    dim=1<<nchan
    X=np.zeros((dim,dim),complex)
    for state in range(dim): X[state^(1<<site),state]=1
    return X


def feshbach(active_ops, multiplicity, t, dm=2.7, b=0.8):
    # multiplicity=1: one channel per bond. multiplicity=2: Kramers closed shell.
    n=active_ops[0].shape[0]
    nchan=len(active_ops)*multiplicity
    cdim=1<<nchan
    I=np.eye(n)
    H0c=np.diag([dm*state.bit_count() for state in range(cdim)])
    V=np.zeros((cdim*n,cdim*n),complex)
    amp=b/math.sqrt(multiplicity)
    ch=0
    for A in active_ops:
        for _ in range(multiplicity):
            V += amp*np.kron(control_x(nchan,ch),A)
            ch+=1
    W4=(b*b/dm)*sum((A@A for A in active_ops),start=np.zeros_like(active_ops[0]))
    H=np.kron(H0c,I)+t*V+t*t*np.kron(np.eye(cdim),W4)
    # control-major ordering; low control is state 0.
    p=np.arange(n)
    q=np.arange(n,cdim*n)
    Hpp=H[np.ix_(p,p)]
    Hpq=H[np.ix_(p,q)]
    Hqq=H[np.ix_(q,q)]
    return herm(Hpp-Hpq@la.solve(Hqq,Hpq.conj().T,assume_a='her'))


def predicted(A,B,multiplicity,dm=2.7,b=0.8):
    comm=A@B-B@A
    local_factor=1.0 if multiplicity==1 else 0.5
    M=local_factor*(A@A@A@A+B@B@B@B)+0.5*comm.conj().T@comm
    return herm((b**4/dm**3)*M)


def main():
    n=4;dm=2.7;b=0.8
    A=herm(RNG.normal(size=(n,n))+1j*RNG.normal(size=(n,n)))
    B=herm(RNG.normal(size=(n,n))+1j*RNG.normal(size=(n,n)))
    A/=opnorm(A);B/=opnorm(B)
    print('KRAMERS CLOSED-SHELL BRIDGE FESHBACH CERTIFICATE')
    print('='*72)
    print('numpy',np.__version__,'scipy',scipy.__version__)
    ok=True
    for mult,name in ((1,'single flavor'),(2,'closed-shell Kramers')):
        target=predicted(A,B,mult,dm,b)
        print('\n'+name)
        prev=None
        for t in (0.18,0.13,0.09,0.065,0.045,0.032,0.025):
            F=feshbach([A,B],mult,t,dm,b)
            ext=F/t**4
            err=opnorm(ext-target)
            print(f'  t={t:.3f}  ||F/t^4-M||={err:.8e}  err/t^2={err/t**2:.8e}')
            if prev is not None: ok &= err<prev
            prev=err
        ok &= prev is not None and prev<2e-6
        print(f'  target norm={opnorm(target):.10e}')
    print('\nThe closed-shell duplicate channel halves only the local A^4+B^4 coefficient; the inter-bond commutator-square coefficient is unchanged.')
    print('OVERALL:','PASS' if ok else 'FAIL')
    return 0 if ok else 1

if __name__=='__main__': raise SystemExit(main())
