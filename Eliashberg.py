#--------------------------------------------------------------------------
# Eliashberg.py - Wrapper file for 2D and 3D eliashberg code allowing c++ class 
# to interface with python scripts
# version: 2
# date: 08/12/2021
# author: Theo Weinberger
#--------------------------------------------------------------------------

from ctypes import *

#import Eliashberg libraries
lib2D = cdll.LoadLibrary('./libeliashberg2D.so')

#set input types for eliashberg functions
lib2D.Eliashberg_New.restype  = c_void_p
lib2D.Eliashberg_New.argtypes = [c_char_p]
lib2D.Eliashberg_SolveEliashberg.restype  = None
lib2D.Eliashberg_SolveEliashberg.argtypes = [c_void_p]

lib3D = cdll.LoadLibrary('./libeliashberg3D.so')

#set input types for eliashberg functions
lib3D.Eliashberg_New.restype  = c_void_p
lib3D.Eliashberg_New.argtypes = [c_char_p]
lib3D.Eliashberg_SolveEliashberg.restype  = None
lib3D.Eliashberg_SolveEliashberg.argtypes = [c_void_p]

#2D Eliashberg Class
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
        self.obj = lib2D.Eliashberg_New(c_config)

    def solve_eliashberg(self):
        """
        Pythonic version of the Eliashberg solver
        """
        lib2D.Eliashberg_SolveEliashberg(self.obj)

#3D Eliashberg Class
class Eliashberg3D(object):
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
        self.obj = lib3D.Eliashberg_New(c_config)

    def solve_eliashberg(self):
        """
        Pythonic version of the Eliashberg solver
        """
        lib3D.Eliashberg_SolveEliashberg(self.obj)
