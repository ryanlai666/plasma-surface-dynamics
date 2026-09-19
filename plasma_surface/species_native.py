"""C++ adapter for the species-resolved conditional sensitivity network."""
import ctypes as ct
from functools import lru_cache
import os,sys
from pathlib import Path
import numpy as np
from .species_kmc import initial_counts,validate
@lru_cache(maxsize=1)
def library():
    suffix='.dll' if os.name=='nt' else '.dylib' if sys.platform=='darwin' else '.so'
    lib=ct.CDLL(str(Path(__file__).resolve().parents[1]/'build'/('species'+suffix)))
    if lib.species_abi_version()!=1:raise RuntimeError('Wrong species ABI')
    f=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');i=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS');l=np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS')
    lib.species_kmc.argtypes=[ct.c_int]*4+[l,i,i,i,f,l,f,ct.c_double,ct.c_uint64,ct.c_uint64,l,l,l];lib.species_kmc.restype=ct.c_int;return lib

def simulate(network,rates,sites,times,seed=0,exposure_end=2.,max_events=10000000):
    validate(network)
    if not isinstance(seed,int) or not 0<=seed<2**64 or not isinstance(max_events,int) or not 0<max_events<2**63:raise ValueError('Invalid seed/event limit')
    rates=np.ascontiguousarray(rates,dtype=np.float64);times=np.ascontiguousarray(times,dtype=np.float64)
    if rates.shape!=(len(network['events']),) or times.ndim!=1 or not len(times):raise ValueError('Invalid array dimensions')
    ns=len(network['states']);ng=len(network['gas_species']);ne=len(rates);gases=list(network['gas_species'])
    src=np.array([e['source'] for e in network['events']],dtype=np.int32);dst=np.array([e['target'] for e in network['events']],dtype=np.int32);ads=np.array([e['kind']=='adsorption' for e in network['events']],dtype=np.int32)
    delta=np.array([[e['gas_delta'].get(g,0) for g in gases] for e in network['events']],dtype=np.int64)
    counts=np.zeros((len(times),ns),dtype=np.int64);gas=np.zeros((len(times),ng),dtype=np.int64);events=np.zeros(ne,dtype=np.int64)
    status=library().species_kmc(ns,ne,ng,len(times),initial_counts(network,sites),src,dst,ads,rates,delta,times,exposure_end,seed,max_events,counts,gas,events)
    if status:raise RuntimeError({1:'Invalid native inputs',2:'Species kMC event budget exceeded',3:'Native runtime failure'}[status])
    return dict(times=times,counts=counts,gas_counts=gas,event_counts=events,events=int(events.sum()),lattice=None)
