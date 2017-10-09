# This code is part of PREFUR

import numpy as np

R = 8.314e-3 # kJ/mol/K

class FES(object):
    """ 
    A class for defining a Free Energy Surface model of protein foldign


    Parameters
    ----------
    protsize : int
        The number of amino acids

    Attributes
    ----------
    natness : array
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

    def __init__(self, protsize):
        self.protsize = protsize
        self.natness = np.linspace(0,1.0,100)
        
    def gen_enthalpy(self, DHres=6.2, kDH=3):
        """
        Generates enthalpy as a function of nativeness

        Parameters
        ----------
        """
        N = self.protsize
        n = self.natness
        self.DHo = N*DHres*(1 + (np.exp(kDH*n) - 1)/(1 - np.exp(kDH)))
        
    def gen_heatcap(self, DCpres=58e-3, kDCp=4.3):
        """
        Generates heat capacity as a function of nativeness

        Parameters
        ----------
        """
        N = self.protsize
        n = self.natness
        self.DCp = N*DCpres*(1 + (np.exp(kDCp*n) - 1)/(1 - np.exp(kDCp)))

    def gen_entropy(self, DSres=16.5e-3):
        """
        Generates entropy as a function of nativeness

        Parameters
        ----------
        DSres : float
            Conformational entropy per residue as a function of nativeness (kJ/mol)

        """
        N = self.protsize
        n = self.natness
        self.DSconf = N*(-R*(n*np.log(n) + (1.-n)*np.log(1-n)) + (1-n)*DSres)
        self.DSconf[0] = N*DSres
        self.DSconf[-1] = 0

    def gen_free(self, T, Tref=385.):
        """
        Generates entropy as a function of nativeness

        Parameters
        ----------
        DSres : float
            Conformational entropy per residue as a function of nativeness (kJ/mol)

        """
        self.DH = self.DHo + self.DCp*(T - Tref)
        self.DS = self.DSconf + self.DCp*np.log(T/Tref)
        self.DG = self.DH - T*self.DS
