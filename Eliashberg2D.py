#--------------------------------------------------------------------------
# Eliashberg2D.py - Wrapper file for 2D eliashberg code allowing c++ class 
# to interface with python scripts
# version: 1
# date: 29/11/2021
# author: Theo Weinberger
#--------------------------------------------------------------------------

from ctypes import *

#import Eliashberg library
lib = cdll.LoadLibrary('./libeliashberg2D.so')

#set input types for eliashberg functions
lib.Eliashberg_New.restype  = c_void_p
lib.Eliashberg_New.argtypes = [c_char_p]
lib.Eliashberg_SolveEliashberg.restype  = None
lib.Eliashberg_SolveEliashberg.argtypes = [c_void_p]

class Eliashberg2D(object):
    """
    Pythonic version of the Eliashberg class
    """

    def __init__(self, config):
        """ Initialiser for the eliashberg class

        Args:

            config: a string containing the name of the Eliashbger configuration file
        """

        #name of configuration file in c code        
        c_config = config.encode('utf-8')

        #create an Eliashberg class object        
        self.obj = lib.Eliashberg_New(c_config)

    def solve_eliashberg(self):
        """
        Pythonic version of the Eliashberg solver
        """
        lib.Eliashberg_SolveEliashberg(self.obj)
