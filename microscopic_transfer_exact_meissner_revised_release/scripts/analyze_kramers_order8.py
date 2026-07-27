#!/usr/bin/env python3
"""Closed-shell Kramers correction to the degree-eight bridge sector.

Two Kramers partners of amplitude b/sqrt(2) reproduce the single-flavor
coefficients through degree six.  At degree eight, partner double excitation
reduces the local sum_e B_e^4 coefficient by one half; the inter-bond
commutator-square coefficient is unchanged.  This script tests whether that
extra term can restore all-filling source safety.
"""
from __future__ import annotations
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as la
import scipy.optimize as opt

ROOT=Path(__file__).parent
spec=importlib.util.spec_from_file_location('tuning',ROOT/'analyze_order8_tuning.py')
tun=importlib.util.module_from_spec(spec);sys.modules[spec.name]=tun;assert spec.loader;spec.loader.exec_module(tun)
mod=tun.mod


def summary(active,S):
    Z=active.Z
    C=mod.herm(Z.conj().T@S)
    X=S-Z@C
    return {"source":mod.opnorm(S),"compression":mod.opnorm(C),"cross":mod.opnorm(X),"eigs":la.eigvalsh(C)}


def serial(x):
    if isinstance(x,np.ndarray): return [float(v) for v in x]
    if isinstance(x,(np.floating,np.integer)): return x.item()
    raise TypeError(type(x).__name__)


def main():
    L,dm,ds,da=3,10.0,8.0,8.0
    s,a=0.9,0.45
    ctrl=mod.control_space(L,dm,ds,da)
    output=[]
    print('KRAMERS CLOSED-SHELL DEGREE-EIGHT AUDIT')
    print('='*78)
    print('The duplicate-partner bridge sector changes')
    print('  V_single = dm^-3 [sum B_e^4 + 1/2 sum_{e<f}[B_e,B_f]^*[B_e,B_f]]')
    print('to')
    print('  V_closed = dm^-3 [1/2 sum B_e^4 + 1/2 sum_{e<f}[B_e,B_f]^*[B_e,B_f]].')
    print()
    for n in (1,2,3):
        active,U,Vsingle,mixed=tun.raw_sources(L,n,ctrl,dm,s,a)
        local=np.zeros_like(active.Z)
        for B in active.B:
            local += (B@B@B@B)@active.Z / dm**3
        Vclosed=Vsingle-0.5*local
        f=lambda x: mod.opnorm(U+x*Vclosed)
        res=opt.minimize_scalar(f,bounds=(0,2),method='bounded',options={'xatol':1e-14})
        x=float(res.x); S=U+x*Vclosed
        data={"pairs":n,"best_x":x,"best_b":math.sqrt(x),"minimum":summary(active,S),"U":summary(active,U),"Vclosed":summary(active,Vclosed),"mixed":mod.opnorm(mixed)}
        output.append(data)
        print(f'n={n}, comps={active.compositions}')
        print(f'  best b                         = {math.sqrt(x):.12e}')
        print(f'  minimum source norm            = {data["minimum"]["source"]:.12e}')
        print(f'  compression norm               = {data["minimum"]["compression"]:.12e}')
        print(f'  cross-source norm              = {data["minimum"]["cross"]:.12e}')
        print('  compression eigenvalues        = ['+', '.join(f'{v:+.8e}' for v in data['minimum']['eigs'])+']')
        print()
    one_exact=output[0]['minimum']['source']<1e-9
    higher_fail=all(x['minimum']['source']>1e-6 for x in output[1:])
    print('one-pair cancellation exists:', 'YES' if one_exact else 'NO')
    print('all-filling Kramers cancellation:', 'NO' if higher_fail else 'UNRESOLVED')
    print('OVERALL:','PASS' if one_exact and higher_fail else 'FAIL')
    (ROOT/'kramers_order8_analysis.json').write_text(json.dumps(output,default=serial,indent=2)+'\n')
    return 0 if one_exact and higher_fail else 1

if __name__=='__main__': raise SystemExit(main())
