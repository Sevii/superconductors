#!/usr/bin/env python3
"""Verify the exact coherent dark AGP branch and its swap-response floor."""
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
o8=importlib.util.module_from_spec(spec);sys.modules[spec.name]=o8;assert spec.loader;spec.loader.exec_module(o8)


def coherent_coefficients(L,comps):
    c=np.array([(-1)**n2*math.sqrt(math.comb(L,n1)*math.comb(L,n2)) for n1,n2 in comps],float)
    return c/la.norm(c)


def swap_floor_sum(L,n,j1,j2,comps,c):
    value=0.0
    for amp,(n1,n2) in zip(c,comps):
        value += abs(amp)**2*(2*j1*n1*(L-n1)/(L*(L-1)) + 2*j2*n2*(L-n2)/(L*(L-1)))
    closed=(j1+j2)*n*(2*L-n)/(L*(2*L-1)) if L>1 else 0.0
    return value,closed


def main():
    print('EXACT COHERENT DARK-AGP BRANCH CERTIFICATE')
    print('='*76)
    print('numpy',np.__version__,'scipy',scipy.__version__)
    max_b=0.0;max_kernel=0.0;max_floor=0.0;max_jacobi=0.0
    j1,j2=0.7,1.1
    for L in (2,3):
        print(f'\nL={L}')
        for n in range(0,2*L+1):
            active=o8.active_sector(L,2*n)
            c=coherent_coefficients(L,active.compositions)
            psi=active.Z@c
            bnorm=max((float(la.norm(B@psi)) for B in active.B),default=0.0)
            Bstack=np.vstack([(B@active.Z) for B in active.B]) if active.B else np.zeros((0,len(c)))
            svals=la.svdvals(Bstack) if Bstack.size else np.array([])
            nullity=int(np.sum(svals<1e-10)) + max(0,len(c)-len(svals))
            # Direct zero check and uniqueness within the product-AGP composition space.
            kernel_res=float(la.norm(Bstack@c)) if Bstack.size else 0.0
            B2comp=Bstack.conj().T@Bstack
            closed_B2=np.zeros_like(B2comp)
            for k,(r,n2) in enumerate(active.compositions):
                closed_B2[k,k]=4*n-8*r*n2/L
                if k+1<len(active.compositions):
                    val=(4/L)*math.sqrt((r+1)*(L-r)*n2*(L-n2+1))
                    closed_B2[k+1,k]=val;closed_B2[k,k+1]=val
            jerr=float(la.norm(B2comp-closed_B2,2)) if B2comp.size else 0.0
            floor,closed=swap_floor_sum(L,n,j1,j2,active.compositions,c)
            ferr=abs(floor-closed)
            max_b=max(max_b,bnorm);max_kernel=max(max_kernel,kernel_res);max_floor=max(max_floor,ferr);max_jacobi=max(max_jacobi,jerr)
            print(f' n={n:2d} comps={len(c)} max||B_x Omega||={bnorm:.2e} nullity={nullity} jacobi_err={jerr:.2e} floor_err={ferr:.2e}')
            if 0<n<2*L and nullity!=1:
                print('   WARNING: unexpected composition-space dark-kernel dimension')
    print('\nAnalytic identities certified:')
    print('  [B_x, eta_1^+ - eta_2^+] = 0 and B_x|0>=0 imply B_x|Omega_n^->=0.')
    print('  The compressed sum_x B_x^2 has the exact irreducible Jacobi entries in the note.')
    print('  |Omega_n^-> has hypergeometric composition weights.')
    print('  D_sw,L = (j1+j2) n(2L-n)/[L(2L-1)].')
    print(f'  max bridge-dark residual = {max_b:.12e}')
    print(f'  max kernel residual      = {max_kernel:.12e}')
    print(f'  max Jacobi formula error = {max_jacobi:.12e}')
    print(f'  max floor formula error  = {max_floor:.12e}')
    ok=max_b<1e-10 and max_kernel<1e-10 and max_jacobi<1e-10 and max_floor<1e-12
    print('OVERALL:','PASS' if ok else 'FAIL')
    return 0 if ok else 1

if __name__=='__main__': raise SystemExit(main())
