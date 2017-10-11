# This code is part of PREFUR

import numpy as np

R = 8.314e-3 # kJ/mol/K

class FES(object):
    """ 
    A class for defining a Free Energy Surface model of protein foldign


    Parameters
    ----------
    nres : int
        The number of amino acids

    Attributes
    ----------
    nat : array
        The nativeness value.

    enthalpy : array
        Enthalpy as a function of nativeness at the reference temperature.

    entropy : array
        Entropy as a function of nativeness at the reference temperature.

    heat_capacity : array
        Heat capacity as a function of nativeness.

    free : array
        Free energy.

    """

    def __init__(self, nres):
        self.nres = nres
        self.nat = np.linspace(0,1.0,100)
        
    def gen_enthalpy(self, DHres=6.2, kDH=3):
        """
        Generates enthalpy as a function of nativeness

        Parameters
        ----------
        DHres : float
            Enthalpy per residue at reference temperature (kJ/mol).

        kDH : float
            Curvature of enthalpy profile.

        """
        N = self.nres
        n = self.nat
        self.DHo = N*DHres*(1 + (np.exp(kDH*n) - 1)/(1 - np.exp(kDH)))
        
    def gen_heatcap(self, DCpres=58e-3, kDCp=4.3):
        """
        Generates heat capacity as a function of nativeness

        Parameters
        ----------
        DCpres : float
            Heat capacity per residue (kJ/mol/K).

        kDCp : float
            Curvature of heat capacity profile.

        """
        N = self.nres
        n = self.nat
        self.DCp = N*DCpres*(1 + (np.exp(kDCp*n) - 1)/(1 - np.exp(kDCp)))

    def gen_entropy(self, DSres=16.5e-3):
        """
        Generates entropy as a function of nativeness

        Parameters
        ----------
        DSres : float
            Conformational entropy per residue as a function of nativeness (kJ/mol)

        """
        N = self.nres
        n = self.nat
        self.DSconf = N*(-R*(n*np.log(n) + (1.-n)*np.log(1-n)) + (1-n)*DSres)
        self.DSconf[0] = N*DSres
        self.DSconf[-1] = 0

    def gen_free(self, temp=298., Tref=385.):
        """
        Generates entropy as a function of nativeness

        Parameters
        ----------
        temp : float
            Temperature of interest (K). 

        Tref : float
            Reference temperature (K). 

        """
        try:
            self.DHo
        except AttributeError:
            self.gen_enthalpy()

        try:
            self.DSconf
        except AttributeError:
            self.gen_entropy()

        try:
            self.DCp
        except AttributeError:
            self.gen_heatcap()

        self.DH = self.DHo + self.DCp*(temp - Tref)
        self.DS = self.DSconf + self.DCp*np.log(temp/Tref)
        self.DG = self.DH - temp*self.DS

    def denature(self, FD, temp=298., Tref=385., C=0.04, j=8.):
        """
        Introduce chemical denaturant

        Parameters
        ----------
        FD : float, array
            Amount of denaturantion free energy.

        Tref : float
            Reference temperature (K). 

        """
        self.mdenat = 1 - (1 + C)*(self.nat**j / (self.nat**j + C))
        self.gen_free(temp, Tref=Tref)

        self.DGdenat = self.DG - self.mdenat*FD

def stability(nat, free, temp=298):
    """
    Calculates the stability given a free energy profile

    Parameters
    ----------
    nat : array
        The nativeness value.

    free : array
        The free energy profile.

    temp : float
        The temperature of interest

    """
    # find transition state
    free_ts = np.max(free[10:-10])
    its = np.argmin(np.abs(free - free_ts))

    # calculate populations
    beta = 1./(temp*8.314e-3) 
    pop = np.exp(-beta*free)
    pop /= np.trapz(pop, nat)
    pf = np.trapz(pop[its:], nat[its:])
    pu = np.trapz(pop[:its], nat[:its])
    
    stab = R*temp*np.log(pf/pu)
    return pf, pu, stab

def barrier(free, n=None, u=None):
    """
    Estimation of barrier height from free energy surface

    Parameters
    ----------
    free : array
        Free energy profile.
        
    """
    fmax = np.max(free[10:-10])
    imax = np.argmin(np.abs(free - fmax))
    fumin = np.min(free[:imax])
    ffmin = np.min(free[imax:])
    return fmax - fumin, fmax - ffmin
