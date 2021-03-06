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
        nres = self.nres
        nat = self.nat
        self.DHo = enthalpy(nres, nat, DHres=DHres, kDH=kDH)

    def gen_enthalpy_local(self, DHres=0, kDH=-1.5):
        """
        Generates local enthalpy as a function of nativeness

        Parameters
        ----------
        DHres : float
            Enthalpy per residue at reference temperature (kJ/mol).

        kDH : float
            Curvature of enthalpy profile.

        """
        nres = self.nres
        nat = self.nat
        self.DHo_loc = enthalpy(nres, nat, DHres=DHres, kDH=kDH)

    def gen_enthalpy_nonlocal(self, DHres=6.2, kDH=3.75):
        """
        Generates non-local enthalpy as a function of nativeness

        Parameters
        ----------
        DHres : float
            Enthalpy per residue at reference temperature (kJ/mol).

        kDH : float
            Curvature of enthalpy profile.

        """
        nres = self.nres
        nat = self.nat
        self.DHo_nonloc = enthalpy(nres, nat, DHres=DHres, kDH=kDH)

    def gen_enthalpy_global(self, DHloc=0, DHnonloc=6.4, kDHloc=-1.5, kDHnonloc=3.75):
        self.gen_enthalpy_local(DHres=DHloc, kDH=kDHloc)
        self.gen_enthalpy_nonlocal(DHres=DHnonloc, kDH=kDHnonloc)
        self.DHo = self.DHo_loc + self.DHo_nonloc
        
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
        nres = self.nres
        nat = self.nat
        self.DCp = nres*DCpres*(1 + (np.exp(kDCp*nat) - 1)/(1 - np.exp(kDCp)))

    def gen_entropy(self, DSres=16.5e-3):
        """
        Generates entropy as a function of nativeness

        Parameters
        ----------
        DSres : float
            Conformational entropy per residue as a function of nativeness (kJ/mol)

        """
        nres = self.nres
        nat = self.nat
        self.DSconf = nres*(-R*(nat*np.log(nat) + (1.-nat)*np.log(1-nat)) + (1-nat)*DSres)
        self.DSconf[0] = nres*DSres
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

def enthalpy(nres, nat, DHres=None, kDH=None):
    """
    Helper function for estimating enthalpies

    Parameters
    ----------
    N : int
        

    n : array
        

    DHres : float
        Enthalpy per residue at reference temperature (kJ/mol).

    kDH : float
        Curvature of enthalpy profile.

    """
    return nres*DHres*(1 + (np.exp(kDH*nat) - 1)/(1 - np.exp(kDH)))

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

def barrier(free):
    """
    Estimation of barrier height from free energy surface

    Parameters
    ----------
    free : array
        Free energy profile.
        
    """
    fmax = np.max(free[30:-10])
    imax = np.argmin(np.abs(free - fmax))
    fumin = np.min(free[:imax])
    ffmin = np.min(free[imax:])
    return fmax - fumin, fmax - ffmin
