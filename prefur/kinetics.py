# This code is part of PREFUR

import numpy as np

def rates(barrier=None, nres=None, temp=298):
    """
    Estimation of (un)folding rates in s^-1 from the free energy barrier

    Parameters
    ----------
    barrier : float
        Free energy barrier in kJ/mol

    nres : int
        The number of amino acids

    T : float
        The temperature of interest
        
    """
    beta = 1/(8.314e-3 * temp)
    return 3.6e6/nres*np.exp(-beta*barrier)
