# This code is part of PREFUR

import numpy as np
from prefur import thermo

def rates(barrier=None, nres=None, temp=298.):
    """
    Estimation of (un)folding rates in s^-1 from the free energy barrier

    Parameters
    ----------
    barrier : float
        Free energy barrier in kJ/mol

    nres : int
        The number of amino acids

    temp : float
        The temperature of interest
        
    """
    beta = 1/(8.314e-3 * temp)
    return 3.6e6/nres*np.exp(-beta*barrier)

def predict(nres=None, struct=None, temp=298.):
    """
    Function for estimating folding and unfolding rates from protein size
    and structural type

    Parameters
    ----------
    nres : int
        The number of amino acids.

    struct : str
        The structural type. Acceptable types are "a", "b" and "ab".

    temp : float
        The temperature of interest.
        
    """
    FES = thermo.FES(nres)
    if not struct:
        pass
    elif struct == "a":
        FES.gen_enthalpy_global(DHloc=2.15, DHnonloc=4.82)
    elif struct == "b":
        FES.gen_enthalpy_global(DHloc=1.31, DHnonloc=5.3)
    elif struct == "ab":
        FES.gen_enthalpy_global(DHloc=1.5, DHnonloc=5.21)
    else:
        print " Unknown structural type "
        print struct
        return
    FES.gen_free(temp=temp) 
    bf, bu = thermo.barrier(FES.DG)
    kf = rates(barrier=bf, nres=nres, temp=temp)
    ku = rates(barrier=bu, nres=nres, temp=temp)
    return kf, ku, FES
