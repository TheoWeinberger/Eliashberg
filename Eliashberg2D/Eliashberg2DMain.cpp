/**
 * @file Eliashberg2DMain.cpp
 * @author Theo Weinberger
 * @brief This code employs the class defined in Eliashberg2D.cpp and Eliashberg.hpp
 * to solve for Tc. It uses the formalism from:
 * 
 *          p-wave and d-wave superconductivity in quasi-two-dimensional metals
 *          P. Monthoux and G. G. Lonzarich
 *          Phys. Rev. B 59, 14598 – Published 1 June 1999
 * 
 * To calculate the energies, frequencies and transition temperatures in MeV and
 * the density of states in 1/meV for a spin mediated superconductor.
 * 
 * It is based of the code Tc_Adaptive.m developed by Ran Tao
 * 
 * This code requires a system to be reduced to the tight-binding approximation 
 * Hamiltonian with dispersion given by
 * 
 *   $\bm{\epsilon_p} = -2t[\cos(p_xa) + \cos(p_ya)] - 4t'\cos(p_xa)\cos(p_ya)
 * 
 * Args:
 * 
 * t: the tight binding hopping matrix element in meV 
 * ratioTight: the ratio of tight-binding hopping elements tPrime/t
 * t0: the initial temperature
 * n0: The initial number of positive fermion Matsubara frequencies within the cutoff
 * omegaC: The cutoff frequency, this should be around 20t but is equal to 2*pi*T0*N(T0)
 * nK: The number of k points along kx, defining the amount of k-point sampling
 * note at this is a 2D system the total number of points goes as Nk^2
 * tSF: The characteristic spin fluctuation temperature. Tsf*kappa02 is approximately constant
 * doping: The doping, calculated by Luttinger's therorem
 * symm: whether the symmetirsation of sigma, the quasiparticle self energy, must be symmetrised
 * chi0: the static susceptibility
 * magModel: Whether the FM or AFM model is being used, takes values FM or AFM
 * phiModel: The model being used for the anomalous self energy, takes values s, p, d
 * errSigma: The relative tolerance for the termination criteria for the convergence in calculating sigma
 * errLambda: The relative tolerance for the termination criteria for the convergence in calculating lambda
 * maxIter: The maximum number of iterations in finding self-consistent solutions at each temperature
 * relaxation: An array containing the mesh of the relaxation weight. This is scanned by the program
 * in search of convergent solutions
 * phiRatio: the ratio between scale of initial trial Phi and dPhi used in solving for lambda
 * plot: Whether we are plotting with respect to kappa or g2chi0t
 * gSquaredchi0tSample: g^2 chi0/t, can be input as an array
 * k0Squarred: k0^2 where k0 is the inverse correlation length without strong magnetic correlations
 * kSquaredSample: where k is the inverse correlation length with strong magnetic correlations, can be an array input
 * 
 * 
 * 
 * 
 * The initial number of postive fermion mastubara frequencies within the cutoff 
 * must all be specified as N_0
 * 
 * 
 * 
 * 
 * 
 * @version 0.1
 * @date 2021-10-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "Eliashberg2D.hpp"
#include "Eliashberg2DSettings.hpp"

/**
 * @brief Main function to call the Eliashberg class and settings
 * 
 * @param argc The number of arguments supplied to the Eliashberg function. This should be 1 as only one configuratations file should be supplied
 * @param argv The arguments supplied to the Eliashberg function, this should be the name of the configurations file
 * @return int 0/1 whether the code has run successfully
 */
int main(int argc, char** argv)
{

    //if more than one input file is supplied, throw an error
    if(argc > 2)
    {

        std::cout << "Too many input files supplied to Eliashberg solver" << std::endl;
        std::cout << "Please only provide one configuration file" << std::endl;
        exit(1);

    }

    //if no input file is supplied, throw an error
    if(argc < 2)
    {

        std::cout << "No input files supplied to Eliashberg solver" << std::endl;
        std::cout << "Please provide a configuration file" << std::endl;
        exit(1);

    }

    //intialise the Eliashberg code with a configurations file
    Eliashberg eliashberg(argv[1]);

    //solve the eliashberg problem
    eliashberg.SolveEliashberg();

    return 0;

}