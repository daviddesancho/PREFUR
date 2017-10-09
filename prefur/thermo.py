# This code is part of PREFUR

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
        
#    def gen_enthalpy(self, DHres=, kDH=kDH):
#        
#    def gen_heatcap(self, DCpres=, kDCp=kDCp):

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
        self.DSconf = [N*(-R*(n*np.log(n) + (1.-n)*np.log(1-n)) + (1-n)*DSres
        self.DSconf[0] = N*DSres
