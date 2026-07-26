#!/usr/bin/env python3
"""Four-band interaction-subbundle screen.

Geometry is evaluated with gauge-invariant projector differences.  Quantities
named Gamma are geometric screening proxies only; any 2*P*Gamma value is a
conditional Jacobi proxy, not a demonstrated physical twist source.
"""
from __future__ import annotations
import math, json, csv, time, argparse, zipfile
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any
import numpy as np

HBAR2_OVER_2ME_MEV_A2 = 3809.98212
A0_A = 3.317

@dataclass(frozen=True)
class ContinuumParams:
    theta_deg: float=3.65
    V_meV: float=9.0
    psi_deg: float=128.0
    w_meV: float=18.0
    mstar: float=0.45
    Vz_meV: float=0.0

@dataclass(frozen=True)
class InteractionProfile:
    name: str
    u_top: float
    u_bottom: float
    inter_x: float
    inter_y: float
    pair_hop_meV: float=2.0

PROFILES=[
 InteractionProfile('symmetric',1.0,1.0,48/35,52/35),
 InteractionProfile('mild_asym',1.0,0.8,48/35,52/35),
 InteractionProfile('orbital_selective',1.0,20/35,48/35,52/35),
]

def hex_indices(shell:int):
    return [(m,n) for m in range(-shell,shell+1) for n in range(-shell,shell+1)
            if max(abs(m),abs(n),abs(m-n))<=shell]

class Model:
    def __init__(self,p:ContinuumParams,shell:int=2):
        self.p=p; self.shell=shell
        self.inds=hex_indices(shell); self.index={x:i for i,x in enumerate(self.inds)}; self.ng=len(self.inds)
        am=A0_A/(2*math.sin(math.radians(p.theta_deg)/2))
        self.aM_A=am
        g=4*math.pi/(math.sqrt(3)*am)
        self.b1=np.array([g,0.]); self.b2=np.array([-g/2,math.sqrt(3)*g/2])
        self.kap_top=(2*self.b1+self.b2)/3
        self.kap_bottom=(self.b1+2*self.b2)/3
        self.G=np.array([m*self.b1+n*self.b2 for m,n in self.inds])
        self.dim=2*self.ng
        self._const=self._build_const()
    def _build_const(self):
        p=self.p; ng=self.ng; H=np.zeros((2*ng,2*ng),complex)
        H[:ng,:ng]+=np.eye(ng)*(p.Vz_meV/2)
        H[ng:,ng:]+=np.eye(ng)*(-p.Vz_meV/2)
        ps=math.radians(p.psi_deg)
        odd=[(1,0),(0,1),(-1,-1)]
        for layer in (0,1):
            s=1 if layer==0 else -1
            off=layer*ng
            for q in odd:
                for qs,amp in ((1,p.V_meV*np.exp(1j*s*ps)),(-1,p.V_meV*np.exp(-1j*s*ps))):
                    dq=(qs*q[0],qs*q[1])
                    for i,(m,n) in enumerate(self.inds):
                        j=self.index.get((m-dq[0],n-dq[1]))
                        if j is not None: H[off+i,off+j]+=amp
        for q in ((0,0),(-1,-1),(0,-1)):
            for i,(m,n) in enumerate(self.inds):
                j=self.index.get((m+q[0],n+q[1]))
                if j is not None:
                    H[i,ng+j]+=p.w_meV; H[ng+j,i]+=p.w_meV
        return H
    def h(self,k):
        H=self._const.copy(); c=HBAR2_OVER_2ME_MEV_A2/self.p.mstar
        for layer,kap in enumerate((self.kap_top,self.kap_bottom)):
            q=k[None,:]+self.G-kap[None,:]
            H[layer*self.ng:(layer+1)*self.ng,layer*self.ng:(layer+1)*self.ng]+=np.diag(-c*np.sum(q*q,axis=1))
        return H
    def eig(self,k,n=4):
        e,u=np.linalg.eigh(self.h(k)); idx=np.argsort(e)[::-1][:n]; return e[idx],u[:,idx]
    def interaction_data(self,k,prof:InteractionProfile,nwin=3):
        e4,u4=self.eig(k,max(4,nwin)); U=u4[:,:nwin]
        A=U[:self.ng,:]; B=U[self.ng:,:]
        Ft=A.conj().T@A; Fb=B.conj().T@B
        C=A.conj().T@B
        Fx=(C+C.conj().T)/math.sqrt(2)
        Fy=(-1j*C+1j*C.conj().T)/math.sqrt(2)
        Fs=(Ft,Fb,Fx,Fy)
        weights=np.array([prof.u_top,prof.u_bottom,prof.inter_x,prof.inter_y],float)
        weights=weights/weights.sum()
        G=np.zeros((nwin,nwin),complex)
        for wt,F in zip(weights,Fs): G+=wt*(F.conj().T@F)
        lam,v=np.linalg.eigh(G); order=np.argsort(lam)[::-1]; lam=lam[order]; v=v[:,order]
        psi=U@v
        hI=v.conj().T@np.diag(e4[:nwin])@v
        srcI=np.zeros((nwin,nwin),float); srcE=np.zeros((nwin,nwin),float)
        for wt,F in zip(weights,Fs):
            FI=v.conj().T@F@v
            srcI+=wt*np.abs(FI)**2
            srcE+=wt*np.abs(F)**2
        return {'e4':e4,'U':U,'lam':lam,'v':v,'psi':psi,'hI':hI,'srcI':srcI,'srcE':srcE}

def kpoints(model,nk=5):
    pts=[]
    for i in range(nk):
      for j in range(nk): pts.append(np.array([(i+0.5)/nk,(j+0.5)/nk]))
    pts += [np.array(x,float) for x in [(0,0),(1/3,1/3),(2/3,1/3),(0.5,0.5),(0,0.5),(0.5,0)]]
    seen=set(); out=[]
    for f in pts:
      key=tuple(np.round(f%1,10))
      if key not in seen: seen.add(key); out.append(f)
    return [(f,f[0]*model.b1+f[1]*model.b2) for f in out]

def norm_cross(S,a,b):
    den=math.sqrt(max(S[a,a]*S[b,b],1e-30)); return float(S[a,b]/den)

def eval_point(p,prof,shell=2,nk=5,geometry=False):
    M=Model(p,shell); datas=[]
    for frac,k in kpoints(M,nk):
        d=M.interaction_data(k,prof); d['frac']=frac; d['k']=k; datas.append(d)
    lams=np.array([d['lam'] for d in datas]); sig=np.sqrt(np.maximum(lams,0))
    h=np.array([d['hI'] for d in datas]); e4=np.array([d['e4'] for d in datas])
    srcI=np.array([d['srcI'] for d in datas]); srcE=np.array([d['srcE'] for d in datas])
    g12=2*(lams[:,0]-lams[:,1])/(lams[:,0]+lams[:,1]+1e-30)
    g23=2*(lams[:,1]-lams[:,2])/(lams[:,1]+lams[:,2]+1e-30)
    leak12=np.array([norm_cross(S,0,1) for S in srcI])
    leak13=np.array([norm_cross(S,0,2) for S in srcI])
    leak23=np.array([norm_cross(S,1,2) for S in srcI])
    eleak12=np.array([norm_cross(S,0,1) for S in srcE])
    eleak13=np.array([norm_cross(S,0,2) for S in srcE])
    eleak23=np.array([norm_cross(S,1,2) for S in srcE])
    diag=np.real(np.diagonal(h,axis1=1,axis2=2))
    bw=np.ptp(diag,axis=0)
    gap34=np.maximum(e4[:,2]-e4[:,3],0)
    ambient_gap=float(np.min(gap34))
    scale=max(ambient_gap,1.0)
    kin12=float(np.sqrt(np.mean(np.abs(h[:,0,1])**2))/scale)
    kina3=float(np.sqrt(np.mean(np.abs(h[:,0,2])**2+np.abs(h[:,1,2])**2))/scale)
    bw_ratio=float(max(bw[0],bw[1])/scale)
    out={**asdict(p),'profile':prof.name,'shell':shell,'nk':nk,'n_kpoints':len(datas),
      'aM_A':M.aM_A,
      'sigma1_mean':float(sig[:,0].mean()),'sigma2_mean':float(sig[:,1].mean()),'sigma3_mean':float(sig[:,2].mean()),
      'sigma1_cv':float(sig[:,0].std()/max(sig[:,0].mean(),1e-30)),
      'sigma2_cv':float(sig[:,1].std()/max(sig[:,1].mean(),1e-30)),
      'gap12_rel_min':float(g12.min()),'gap12_rel_mean':float(g12.mean()),
      'gap23_rel_min':float(g23.min()),'gap23_rel_mean':float(g23.mean()),
      'source_leak12_mean':float(leak12.mean()),'source_leakA3_mean':float(np.sqrt((leak13**2+leak23**2)/2).mean()),
      'energy_source_leak12_mean':float(eleak12.mean()),'energy_source_leakA3_mean':float(np.sqrt((eleak13**2+eleak23**2)/2).mean()),
      'kinetic_mix12_over_gap34':kin12,'kinetic_mixA3_over_gap34':kina3,
      'block1_diag_bandwidth_meV':float(bw[0]),'block2_diag_bandwidth_meV':float(bw[1]),'remote_diag_bandwidth_meV':float(bw[2]),
      'ambient_gap34_min_meV':ambient_gap,'active_bandwidth_over_gap34':bw_ratio,
      'energy_top_bandwidth_meV':float(np.ptp(e4[:,0])),'energy_second_bandwidth_meV':float(np.ptp(e4[:,1])),'energy_third_bandwidth_meV':float(np.ptp(e4[:,2])),
    }
    if geometry:
        hstep=1e-3*np.linalg.norm(M.b1)
        metrics=np.zeros((3,2)); gamma=np.zeros((2,2)); count=0
        for frac,k in kpoints(M,max(4,nk)):
            d0=M.interaction_data(k,prof); P=[np.outer(d0['psi'][:,a],d0['psi'][:,a].conj()) for a in range(3)]
            for i in (0,1):
                sh=np.zeros(2); sh[i]=hstep
                dp=M.interaction_data(k+sh,prof); dm=M.interaction_data(k-sh,prof)
                for a in range(3):
                    Pp=np.outer(dp['psi'][:,a],dp['psi'][:,a].conj()); Pm=np.outer(dm['psi'][:,a],dm['psi'][:,a].conj())
                    dP=(Pp-Pm)/(2*hstep)
                    metrics[a,i]+=0.5*np.trace(dP@dP).real
                    if a<2: gamma[a,i]+=np.trace(P[2]@dP@dP).real
            count+=1
        metrics/=count; gamma/=count
        metrics_dimless=metrics/(M.aM_A**2); gamma_dimless=gamma/(M.aM_A**2)
        out.update({
          'metric_block1_x_aM2':float(metrics_dimless[0,0]),'metric_block1_y_aM2':float(metrics_dimless[0,1]),
          'metric_block2_x_aM2':float(metrics_dimless[1,0]),'metric_block2_y_aM2':float(metrics_dimless[1,1]),
          'gamma_1to3_x_aM2':float(gamma_dimless[0,0]),'gamma_1to3_y_aM2':float(gamma_dimless[0,1]),
          'gamma_2to3_x_aM2':float(gamma_dimless[1,0]),'gamma_2to3_y_aM2':float(gamma_dimless[1,1]),
        })
        gx=float((gamma_dimless[0,0]+gamma_dimless[1,0])/2); gy=float((gamma_dimless[0,1]+gamma_dimless[1,1])/2)
        out['Gamma_x_proxy_aM2']=gx; out['Gamma_y_proxy_aM2']=gy
        out['Gamma_x_aM2']=gx; out['Gamma_y_aM2']=gy
        out['conditional_t0_proxy_x_meV_aM2']=2*gx*prof.pair_hop_meV
        out['conditional_t0_proxy_y_meV_aM2']=2*gy*prof.pair_hop_meV
        out['one_pair_hopping_proxy_x_meV_aM2']=2*gx*prof.pair_hop_meV
        out['one_pair_hopping_proxy_y_meV_aM2']=2*gy*prof.pair_hop_meV
    cv=out['sigma1_cv']+out['sigma2_cv']; leak=out['source_leak12_mean']+out['source_leakA3_mean']
    # Log-compressed control penalties; score is only a ranking diagnostic.
    score=(1.4*out['gap12_rel_min']+1.8*out['gap23_rel_min']-1.2*cv-1.1*leak
           -0.8*math.log1p(kin12)-0.9*math.log1p(kina3)-0.45*math.log1p(bw_ratio)
           +0.20*math.log1p(max(ambient_gap,0)/5))
    out['screen_score']=float(score)
    out['qgn_error_max']=float(max(out['sigma1_cv'],out['sigma2_cv'],out['source_leak12_mean'],out['source_leakA3_mean']))
    out['control_error_max']=float(max(kin12,kina3,bw_ratio))
    return out

def interaction_data_window(model: Model, k, prof: InteractionProfile, nwin: int=4):
    eall,uall=model.eig(k,nwin+1); U=uall[:,:nwin]
    A=U[:model.ng,:]; B=U[model.ng:,:]
    Ft=A.conj().T@A; Fb=B.conj().T@B
    C=A.conj().T@B
    Fx=(C+C.conj().T)/math.sqrt(2)
    Fy=(-1j*C+1j*C.conj().T)/math.sqrt(2)
    Fs=(Ft,Fb,Fx,Fy)
    weights=np.array([prof.u_top,prof.u_bottom,prof.inter_x,prof.inter_y],float)
    weights=weights/weights.sum()
    G=np.zeros((nwin,nwin),complex)
    for wt,F in zip(weights,Fs): G+=wt*(F.conj().T@F)
    lam,v=np.linalg.eigh(G); order=np.argsort(lam)[::-1]; lam=lam[order]; v=v[:,order]
    psi=U@v
    hI=v.conj().T@np.diag(eall[:nwin])@v
    srcI=np.zeros((nwin,nwin),float); srcE=np.zeros((nwin,nwin),float)
    for wt,F in zip(weights,Fs):
        FI=v.conj().T@F@v
        srcI+=wt*np.abs(FI)**2
        srcE+=wt*np.abs(F)**2
    return {'eall':eall,'U':U,'lam':lam,'v':v,'psi':psi,'hI':hI,'srcI':srcI,'srcE':srcE}

def eval_point_window(p,prof,shell=2,nk=5,nwin=4,nactive=2,geometry=False):
    if nwin <= nactive:
        raise ValueError("nwin must exceed nactive")
    M=Model(p,shell); datas=[]
    for frac,k in kpoints(M,nk):
        d=interaction_data_window(M,k,prof,nwin); d['frac']=frac; d['k']=k; datas.append(d)
    lams=np.array([d['lam'] for d in datas]); sig=np.sqrt(np.maximum(lams,0))
    h=np.array([d['hI'] for d in datas]); eall=np.array([d['eall'] for d in datas])
    srcI=np.array([d['srcI'] for d in datas]); srcE=np.array([d['srcE'] for d in datas])
    g12=2*(lams[:,0]-lams[:,1])/(lams[:,0]+lams[:,1]+1e-30)
    g2r=2*(lams[:,nactive-1]-lams[:,nactive])/(lams[:,nactive-1]+lams[:,nactive]+1e-30)
    def cross_arrays(Sarr):
        l12=np.array([norm_cross(S,0,1) for S in Sarr])
        lar=[]
        for S in Sarr:
            vals=[norm_cross(S,a,r) for a in range(nactive) for r in range(nactive,nwin)]
            lar.append(math.sqrt(float(np.mean(np.square(vals)))))
        return l12,np.array(lar)
    leak12,leakAR=cross_arrays(srcI); eleak12,eleakAR=cross_arrays(srcE)
    diag=np.real(np.diagonal(h,axis1=1,axis2=2)); bw=np.ptp(diag,axis=0)
    gapout=np.maximum(eall[:,nwin-1]-eall[:,nwin],0)
    ambient_gap=float(np.min(gapout)); scale=max(ambient_gap,1.0)
    kin12=float(np.sqrt(np.mean(np.abs(h[:,0,1])**2))/scale)
    kinAR=float(np.sqrt(np.mean(np.sum(np.abs(h[:,:nactive,nactive:])**2,axis=(1,2))))/scale)
    bw_ratio=float(max(bw[:nactive])/scale)
    out={**asdict(p),'profile':prof.name,'shell':shell,'nk':nk,'nwin':nwin,'nactive':nactive,
      'n_kpoints':len(datas),'aM_A':M.aM_A,
      'sigma1_mean':float(sig[:,0].mean()),'sigma2_mean':float(sig[:,1].mean()),
      'sigma_remote1_mean':float(sig[:,nactive].mean()),
      'sigma1_cv':float(sig[:,0].std()/max(sig[:,0].mean(),1e-30)),
      'sigma2_cv':float(sig[:,1].std()/max(sig[:,1].mean(),1e-30)),
      'gap12_rel_min':float(g12.min()),'gap12_rel_mean':float(g12.mean()),
      'gap2R_rel_min':float(g2r.min()),'gap2R_rel_mean':float(g2r.mean()),
      'source_leak12_mean':float(leak12.mean()),'source_leakAR_mean':float(leakAR.mean()),
      'energy_source_leak12_mean':float(eleak12.mean()),'energy_source_leakAR_mean':float(eleakAR.mean()),
      'kinetic_mix12_over_gapout':kin12,'kinetic_mixAR_over_gapout':kinAR,
      'block1_diag_bandwidth_meV':float(bw[0]),'block2_diag_bandwidth_meV':float(bw[1]),
      'ambient_gap_out_min_meV':ambient_gap,'active_bandwidth_over_gapout':bw_ratio,
      'energy_top_bandwidth_meV':float(np.ptp(eall[:,0])),'energy_second_bandwidth_meV':float(np.ptp(eall[:,1])),
    }
    if geometry:
        hstep=1e-3*np.linalg.norm(M.b1)
        metrics=np.zeros((nactive,2)); gamma=np.zeros((nactive,2)); count=0
        for frac,k in kpoints(M,max(4,nk)):
            d0=interaction_data_window(M,k,prof,nwin)
            P=[np.outer(d0['psi'][:,a],d0['psi'][:,a].conj()) for a in range(nwin)]
            PR=sum(P[nactive:])
            for i in (0,1):
                sh=np.zeros(2); sh[i]=hstep
                dp=interaction_data_window(M,k+sh,prof,nwin)
                dm=interaction_data_window(M,k-sh,prof,nwin)
                for a in range(nactive):
                    Pp=np.outer(dp['psi'][:,a],dp['psi'][:,a].conj())
                    Pm=np.outer(dm['psi'][:,a],dm['psi'][:,a].conj())
                    dP=(Pp-Pm)/(2*hstep)
                    metrics[a,i]+=0.5*np.trace(dP@dP).real
                    gamma[a,i]+=np.trace(PR@dP@dP).real
            count+=1
        metrics/=count; gamma/=count
        md=metrics/(M.aM_A**2); gd=gamma/(M.aM_A**2)
        out.update({
          'metric_block1_x_aM2':float(md[0,0]),'metric_block1_y_aM2':float(md[0,1]),
          'metric_block2_x_aM2':float(md[1,0]),'metric_block2_y_aM2':float(md[1,1]),
          'gamma_1toR_x_aM2':float(gd[0,0]),'gamma_1toR_y_aM2':float(gd[0,1]),
          'gamma_2toR_x_aM2':float(gd[1,0]),'gamma_2toR_y_aM2':float(gd[1,1]),
        })
        gx=float((gd[0,0]+gd[1,0])/2); gy=float((gd[0,1]+gd[1,1])/2)
        out['Gamma_x_proxy_aM2']=gx; out['Gamma_y_proxy_aM2']=gy
        out['Gamma_x_aM2']=gx; out['Gamma_y_aM2']=gy
        out['conditional_t0_proxy_x_meV_aM2']=2*gx*prof.pair_hop_meV
        out['conditional_t0_proxy_y_meV_aM2']=2*gy*prof.pair_hop_meV
        out['one_pair_hopping_proxy_x_meV_aM2']=2*gx*prof.pair_hop_meV
        out['one_pair_hopping_proxy_y_meV_aM2']=2*gy*prof.pair_hop_meV
    cv=out['sigma1_cv']+out['sigma2_cv']; leak=out['source_leak12_mean']+out['source_leakAR_mean']
    score=(1.4*out['gap12_rel_min']+1.8*out['gap2R_rel_min']-1.2*cv-1.1*leak
           -0.8*math.log1p(kin12)-0.9*math.log1p(kinAR)-0.45*math.log1p(bw_ratio)
           +0.20*math.log1p(max(ambient_gap,0)/5))
    out['screen_score']=float(score)
    out['qgn_error_max']=float(max(out['sigma1_cv'],out['sigma2_cv'],out['source_leak12_mean'],out['source_leakAR_mean']))
    out['control_error_max']=float(max(kin12,kinAR,bw_ratio))
    return out

# -----------------------------------------------------------------------------
# Reproducible command-line driver
# -----------------------------------------------------------------------------

VALIDATED_RANGES = {
    'theta_deg': (2.0, 5.0),
    'V_meV': (6.0, 12.0),
    'psi_deg': (112.0, 142.0),
    'w_meV': (12.0, 24.0),
    'mstar': (0.38, 0.55),
    'Vz_meV': (0.0, 25.0),
}

EXPLORATORY_RANGES = {
    'theta_deg': (1.2, 3.2),
    'V_meV': (4.0, 15.0),
    'psi_deg': (100.0, 150.0),
    'w_meV': (10.0, 32.0),
    'mstar': (0.30, 0.70),
    'Vz_meV': (0.0, 40.0),
}

REFERENCE_CASES = {
    'published_3p65': ContinuumParams(3.65, 9.0, 128.0, 18.0, 0.45, 0.0),
    'best_validated_2to5': ContinuumParams(
        2.0249742617785738, 8.942875735250286, 124.17739453787478,
        23.25844058808712, 0.5035757685318687, 3.4914346067966604),
    'balanced_exploratory': ContinuumParams(
        1.2835427677189872, 8.800868865081524, 102.85138920929685,
        28.535289241798942, 0.5714708424613053, 23.317285681012464),
    'deep_control_small_gap': ContinuumParams(
        1.2053778262930535, 14.727734931059292, 127.0012352190945,
        18.37307309597066, 0.6514710306404935, 33.83527386987681),
}


def _lhs_params(n: int, ranges: dict[str, tuple[float, float]], seed: int) -> list[ContinuumParams]:
    """Generate a deterministic Latin-hypercube parameter list."""
    try:
        from scipy.stats import qmc
        sampler = qmc.LatinHypercube(d=len(ranges), seed=seed)
        unit = sampler.random(n)
    except Exception:
        rng = np.random.default_rng(seed)
        unit = np.empty((n, len(ranges)))
        for j in range(len(ranges)):
            bins = (np.arange(n) + rng.random(n)) / n
            rng.shuffle(bins)
            unit[:, j] = bins
    keys = list(ranges)
    lo = np.array([ranges[k][0] for k in keys])
    hi = np.array([ranges[k][1] for k in keys])
    values = lo + unit * (hi - lo)
    return [ContinuumParams(**{k: float(v) for k, v in zip(keys, row)}) for row in values]


def _structured_reference_grid() -> list[ContinuumParams]:
    """Angle/displacement-field grid around the conventional 3.65-degree model."""
    angles = sorted(set(np.round(np.linspace(2.0, 5.0, 31), 10).tolist() + [3.65]))
    fields = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0]
    return [ContinuumParams(theta, 9.0, 128.0, 18.0, 0.45, vz)
            for theta in angles for vz in fields]


def _write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    if not rows:
        raise ValueError('No rows to write')
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({key for row in rows for key in row})
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _evaluate_many(
    specs: list[tuple[ContinuumParams, InteractionProfile]],
    *, shell: int, nk: int, geometry: bool, jobs: int,
) -> list[dict[str, Any]]:
    def work(spec: tuple[ContinuumParams, InteractionProfile]) -> dict[str, Any]:
        p, profile = spec
        return eval_point_window(p, profile, shell=shell, nk=nk, nwin=4, nactive=2,
                                 geometry=geometry)
    if jobs == 1:
        return [work(spec) for spec in specs]
    try:
        from joblib import Parallel, delayed
        return Parallel(n_jobs=jobs, verbose=5)(delayed(work)(spec) for spec in specs)
    except Exception as exc:
        print(f'Parallel execution unavailable ({exc}); falling back to serial.')
        return [work(spec) for spec in specs]


def run_scan(
    *, ranges: dict[str, tuple[float, float]], n_lhs: int, seed: int,
    output_csv: Path, shell: int=2, nk: int=4, jobs: int=1,
    include_structured: bool=False,
) -> list[dict[str, Any]]:
    params = _lhs_params(n_lhs, ranges, seed)
    if include_structured:
        params.extend(_structured_reference_grid())
    specs = [(p, profile) for p in params for profile in PROFILES]
    rows = _evaluate_many(specs, shell=shell, nk=nk, geometry=False, jobs=jobs)
    _write_csv(rows, output_csv)
    return rows


def reproduce_reference_cases(output_json: Path, shell: int=4, nk: int=12) -> dict[str, Any]:
    profile = next(p for p in PROFILES if p.name == 'symmetric')
    results = {
        name: eval_point_window(p, profile, shell=shell, nk=nk, nwin=4, nactive=2,
                                geometry=True)
        for name, p in REFERENCE_CASES.items()
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(results, indent=2))
    return results


def _main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            'Scan a first-harmonic tWSe2 continuum family using the two leading '
            'interaction singular subbundles as QGN block candidates.'))
    parser.add_argument(
        '--mode', choices=['quick', 'validated', 'exploratory', 'reference'],
        default='reference')
    parser.add_argument('--output', type=Path,
                        default=Path('wse2_interaction_subbundle_output'))
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--seed', type=int, default=20260725)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    if args.mode == 'reference':
        path = args.output / 'reference_cases.json'
        reproduce_reference_cases(path)
        print(path)
        return

    if args.mode == 'quick':
        run_scan(ranges=VALIDATED_RANGES, n_lhs=30, seed=args.seed,
                 output_csv=args.output/'quick_scan.csv', shell=2, nk=3,
                 jobs=args.jobs, include_structured=False)
    elif args.mode == 'validated':
        run_scan(ranges=VALIDATED_RANGES, n_lhs=900, seed=args.seed,
                 output_csv=args.output/'validated_scan.csv', shell=2, nk=4,
                 jobs=args.jobs, include_structured=True)
    else:
        run_scan(ranges=EXPLORATORY_RANGES, n_lhs=1800, seed=args.seed+1,
                 output_csv=args.output/'exploratory_scan.csv', shell=2, nk=4,
                 jobs=args.jobs, include_structured=False)


if __name__ == '__main__':
    _main()
