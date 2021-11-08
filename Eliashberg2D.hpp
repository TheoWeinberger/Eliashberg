/**
 * @file Eliashberg2D.hpp
 * @author Theo Weinberger
 * @brief This class uses the formalism from:
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

#include <iostream>
#include <armadillo>
#include <cmath>
#include <complex.h>
#include <fftw3.h>
#include <libconfig.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <stdio.h>
#include <math.h>
#include <libalglib/interpolation.h>

#ifndef ELIASHBERG2D_H
#define ELIASHBERG2D_H


/**
 * @brief A class object that is used to solve the Eliashberg equations for given 
 * inputs
 * 
 */
class Eliashberg
{
public:


    /**
     * @brief Construct a new Eliashberg object using the default paramters
     * 
     */
    Eliashberg();

    /**
     * @brief Construct a new Eliashberg object from input data from configuration 
     * file called, eliashberg.cfg
     * 
     */
    Eliashberg(const std::string&);

    /**
     * @brief Solves the system according to the formalsim in P. Monthoux and G. G. Lonzarich (1999)
     * 
     */
    void SolveEliashberg();

    /**
     * @brief Get the value of tC for the system
     * 
     * @return tC: The critical temperature for superconductivity
     */
    double GetTC()const;


private:

    /**
     * @brief Function that calculates the energy of the system from its
     * dispersion relationsip
     * 
     * @param pX Momentum space vectos in the x direction
     * @param pY Momentum space vectos in the y direction
     * @return dispersion: a vector containing the dispersion relation of the system as a function of qX and qY
     */
    arma::mat _Dispersion(const arma::vec&, const arma::vec&);

    /**
     * @brief Function that calculates the negative root of momentum 
     * part of the susceptibility
     * 
     * @param qX Momentum space vectos in the x direction
     * @param qY Momentum space vectos in the y direction
     * @return qHatMinus: a vector containing the negative root of the susceptibility pole
     */
    arma::vec _QMinus(const arma::vec&, const arma::vec&);

    /**
     * @brief Function that calculates the positive root of momentum 
     * part of the susceptibility
     * 
     * @param qX Momentum space vectos in the x direction
     * @param qY Momentum space vectos in the y direction
     * @return qHatPlus: a vector containing the positive root of the susceptibility pole
     */
    arma::vec _QPlus(const arma::vec&, const arma::vec&);

    /**
     * @brief Function that calculates the eta parameter for the suscpetibility
     * eta is a momentum dependent coherence length (I think)
     * 
     * @param qX Momentum space vectos in the x direction
     * @param qY Momentum space vectos in the y direction
     */
    void _Eta(const arma::vec&, const arma::vec&);

    /**
     * @brief Function that calculates the cutoff frequency for the dynamics susceptibility integration
     * 
     * @return omega0: The cutoff frequency for the integration of the dynamic susceptibility
     */
    arma::vec _Omega0();

        /**
     * @brief Function that calculates the negative root of momentum 
     * part of the susceptibility. This version is for iterative calculations
     * 
     * @param qX Momentum space vectos in the x direction
     * @param qY Momentum space vectos in the y direction
     * @return qHatMinus: a vector containing the negative root of the susceptibility pole
     */
    double _QMinus(const double&, const double&);

    /**
     * @brief Function that calculates the positive root of momentum 
     * part of the susceptibility. This version is for iterative calculations
     * 
     * @param qX Momentum space vectos in the x direction
     * @param qY Momentum space vectos in the y direction
     * @return qHatPlus: a vector containing the positive root of the susceptibility pole
     */
    double _QPlus(const double&, const double&);

    /**
     * @brief Function that calculates the eta parameter for the suscpetibility
     * eta is a momentum dependent coherence length (I think)/. 
     * This version is for iterative calculations
     * 
     * @param qX Momentum space vectos in the x direction
     * @param qY Momentum space vectos in the y direction
     */
    double _Eta(const double&, const double&);

    /**
     * @brief Function that calculates the cutoff frequency for the dynamics susceptibility integration.
     * This version is for iterative calculations
     * 
     * @param qX Momentum space vectors in the x direction
     * @param qY Momentum space vectors in the y direction
     * @return omega0: The cutoff frequency for the integration of the dynamic susceptibility
     */
    double _Omega0(const double&, const double&);

    /**
     * @brief Function to solve for the number of fermions in the system
     * 
     * @param mu the chemical potential
     * @param t the temperature
     * @return muSolve: the solved for value of mu
     */
    arma::mat _nFermi(const double&, const double&);

    /**
     * @brief Function used to evaluate the total state density
     * 
     * @param mu the chemical potential
     * @param p parameters for the equations - the temperature
     * @return muSolve: the solved for value of mu
     */
    double _nTotalEval(double, void*);

    /**
     * @brief Function to solve an eigenvalue equation to find potential 
     * values of mu
     * 
     * @param tInit Initial temperature
     * @return muVec: A vector containing the values of mu
     */
    arma::vec _calcMu(const double&);

    /**
     * @brief Function to calculate the integral of the dynamic susceptibility
     * 
     * @param vmatsu The vector of matusbara frequencies
     * @param qx The mopmentum space vector in the x direction
     * @param qy The momentum space vector in the y direction
     * @return chiQNu The integrated dynamics susceptibility
     */
    arma::cube _ChiQNu(const arma::vec&, const arma::vec&, const arma::vec&);

    /**
     * @brief Function for the analytical expression of the dynamic susceptibility integral
     * 
     * @param y the variable that is being integrated over
     * @param a a system parameter
     * @param b a system parameter
     * @return The value of the integral at y
     */
    double _ChiInt(const double&, const double&, const double&);

    /**
     * @brief This function uses the FFT to calculate a circular convolution between matrix A and B in 
     * k-space. This is then padded with zeros along the highest dimesnion to get the linear convolution
     * in frequency. The inputs are thre dimension arrays and the rnage of fequencies and wavevectors
     * are cetnreed around zero. The output is truncated at the highest dimension from the low index
     * to the high index
     * 
     * @param matrixA The first input matrix
     * @param matrixB The second input matrix
     * @param lowindex The low index for the truncation
     * @param highIndex The high index for the truncation
     * @param firstRun boolean specifying if it is the first run
     * @param rg boolean specifying whther its rg corrections
     * @return matrixConv The matrix corresponding to a linear convolution in frequency of A and B
     */
    arma::cx_cube _MatsuConv(const arma::cx_cube&, const arma::cx_cube&, const int&, const int&, const bool&, const bool&);

    /**
     * @brief Function to set DFT plans for the matsubara frequency convolution
     * 
     * @param in The matrix being transformed
     * @param out The output matrix of the DFT
     */
    void _SetDFTPlans(const arma::cx_cube&, const arma::cx_cube&);

    /**
     * @brief Function to set DFT plans for the matsubara frequency convolution
     * 
     * @param in The matrix being transformed
     * @param out The output matrix of the DFT
     */
    void _SetDFTPlansRG(const arma::cx_cube& in, const arma::cx_cube& out);

    /**
     * @brief Function to delete DFT plans for the matsubara frequency convolution
     * 
     */
    void _DeleteDFTPlans();

    /**
     * @brief Calculate the phi function used to iteratively calculate the eigenalue lambda
     * via the power method
     * 
     * @param qX A vector containing the momentum in the x direction
     * @param qY A vector conaining the momentum in the y direction
     * @return phi A cube containing the data for phi
     */
    arma::cube _PhiFun(const arma::vec&, const arma::vec&);

    /**
     * @brief Calculate the symmetric phi function used to iteratively calculate the eigenalue lambda
     * via the power method
     * 
     * @param qX A vector containing the momentum in the x direction
     * @param qY A vector conaining the momentum in the y direction
     * @return phi A cube containing the data for symmetric phi
     */
    arma::cube _PhiSymm(const arma::vec&, const arma::vec&);

    /**
     * @brief A function to symmetrise a matrix according to its symmetry model (s, p ,d)
     * and the sign symmetry of the matrix filter. The size of the matrix should be even
     * in the wavevector dimensions and so the filer is calculated via the _PhiSymm function.
     * The symmetries are the same as in _PhiFun.
     * 
     * @param matrixA The input matrix
     * @param filter The symmetry filter
     * @return matrix S The output symmetrised matrix
     */
    arma::cube _SymmByFiltLabel(arma::cube&, const arma::cube&);

    /**
     * @brief t: the tight binding hopping matrix element in meV 
     * 
     */
    double _t;

    /**
     * @brief ratioTight: the ratio of tight-binding hopping elements tPrime/t
     * 
     */
    double _ratioTight;

    /**
     * @brief tPrime: tight binding hopping matrix elements in meV
     * 
     */
    double _tPrime;

    /**
     * @brief a: The side length of a square lattice unit cell
     * 
     */
    double _a;

    /**
     * @brief t0: Initial temperature
     * 
     */
    double _t0;

    /**
     * @brief n0: The initial number of positive fermion Matsubara frequencies within the cutoff
     * 
     */
    int _n0;

    /**
     * @brief omegaC: The cutoff frequency, this should be around 20t but is equal to 2*pi*T0*N(T0)
     * 
     */
    double _omegaC;

    /**
     * @brief nK: The number of k points along kx, defining the amount of k-point sampling
     * note at this is a 2D system the total number of points goes as Nk^2
     * 
     */
    int _nK;

    /**
     * @brief tSF: The characteristic spin fluctuation temperature. Tsf*kappa02 is approximately constant
     * 
     */
    double _tSF;

    /**
     * @brief doping: The doping, calculated by Luttinger's therorem
     * 
     */
    double _doping;

    /**
     * @brief symm: whether the symmetirsation of sigma, the quasiparticle self energy, must be symmetrised
     * 
     */
    bool _symm;

    /**
     * @brief chi0: the static susceptibility
     * 
     */
    double _chi0;

    /**
     * @brief magModel: Whether the FM or AFM model is being used, takes values FM or AFM
     * 
     */
    std::string _magModel;

    /**
     * @brief phiModel:The model being used for the anomalous self energy, takes values s, p, d
     * 
     */
    std::string _phiModel;

    /**
     * @brief errSigma: The relative tolerance for the termination criteria for the convergence in calculating sigma
     * 
     */
    double _errSigma;

    /**
     * @brief errLambda: The relative tolerance for the termination criteria for the convergence in calculating lambda
     * 
     */
    double _errLambda;

    /**
     * @brief maxIter: The maximum number of iterations in finding self-consistent solutions at each temperature
     * 
     */
    int _maxIter;

    /**
     * @brief relaxation: An array containing the mesh of the relaxation weight. This is scanned by the program
     * in search of convergent solutions
     * 
     */
    arma::vec _relaxation;

    /**
     * @brief phiRatio: the ratio between scale of initial trial Phi and dPhi used in solving for lambda
     * 
     */
    double _phiRatio;

    /**
     * @brief plot: Whether we are plotting with respect to k or g2chi0t
     * takes values of k or g
     * 
     */
    std::string _plot;

    /**
     * @brief gSquaredchi0tSample: g^2 chi0/t, can be input as an array
     * 
     */
    arma::vec _gSquaredChi0tSample;

    /**
     * @brief k0Squared: k0^2 where k0 is the inverse correlation length without strong magnetic correlations
     * 
     */
    double _k0Squared;

    /**
     * @brief kSquaredSample: where k is the inverse correlation length with strong magnetic correlations, can be an array input
     * 
     */
    arma::vec _kSquaredSample;

    /**
     * @brief kSquared: The actual value of kSquared, talen from kSquaredSample
     * 
     */
    double _kSquared;

    /**
     * @brief Eta: a momentum dependent coherence length (I think)
     * 
     */
    arma::vec _eta;

    /**
     * @brief The q space vector being used, depends on whether FM or AFM
     * 
     */
    arma::vec _qHat;

    /**
     * @brief The chemical potential
     * 
     */
    double _mu;

    /**
     * @brief The 2D dispersion of the system
     * 
     */
    arma::mat _energy;

    /**
     * @brief Max number of iterations when solving for mu
     * 
     */
    int _maxMuIter;

    /**
     * @brief The plan for forward FFTs in the matsubara convolution function
     * 
     */
    fftw_plan _forwardPlan;

    /**
     * @brief The plan for inverse FFTs in the matsubara convolution function
     * 
     */
    fftw_plan _inversePlan;

    /**
     * @brief The plan for forward FFTs in the matsubara convolution function for RG corrections
     * 
     */
    fftw_plan _forwardPlanRG;

    /**
     * @brief The plan for inverse FFTs in the matsubara convolution function for RG corrections
     * 
     */
    fftw_plan _inversePlanRG;



};


/**
 * @brief Function to solve for the number of fermions in the system
 * 
 * @param mu the chemical potential
 * @param t the temperature
 * @return arma::mat nFermiMat matrix of the nFermi at each sampling point
 */
arma::mat NFermi(const double& mu, const double& t, const arma::mat& energy)
{

    arma::mat nFermiMat = 1.0/(exp((energy-mu)/t) + 1);

    return nFermiMat;

}

/**
 * @brief Function to solve for the derivate of the
 * number of fermions in the system
 * 
 * @param mu the chemical potential
 * @param t the temperature
 * @return double nFermiMat matrix of the derivative of nFermi at each sampling point
 */
arma::mat NFermiDeriv(const double& mu, const double& t, const arma::mat& energy)
{

    arma::mat nFermiMatDeriv = exp((energy-mu)/t)/(t*pow((exp((energy-mu)/t) + 1), 2.0));

    return nFermiMatDeriv;

}

/**
 * @brief Parameter structure for solving nTotal
 * 
 */
struct NTotalEvalParams
{ 
    //temperature
    double t; 

    //energy
    arma::mat energy; 

    //k space sampling
    int nK; 

    //doping by Luttinger's theorem
    double doping;
};

/**
 * @brief Function used to evaluate the total state density
 * 
 * @param mu the chemical potential
 * @param p parameters for the equations - the temperature
 * @return nTotal the total density of states
 */
double NTotal(double mu, void* p)
{

    struct NTotalEvalParams * params = (struct NTotalEvalParams *)p;

    double t = (params->t);
    arma::mat energy = (params->energy);
    double nK = (params->nK); 
    double doping = (params->doping);

    double nTotal = 2.0*arma::accu(NFermi(mu, t, energy))/pow(nK, 2.0) - doping;

    return nTotal;
}

/**
 * @brief Function used to evaluate the derivative of the total state density
 * 
 * @param mu the chemical potential
 * @param p parameters for the equations - the temperature
 * @return nTotalDeriv the derivative of the total density of states
 */
double NTotalDeriv(double mu, void* p)
{

    struct NTotalEvalParams * params = (struct NTotalEvalParams *)p;

    double t = (params->t);
    arma::mat energy = (params->energy);
    double nK = (params->nK); 

    double nTotalDeriv = 2.0*arma::accu(NFermiDeriv(mu, t, energy))/pow(nK, 2.0);

    return nTotalDeriv;
}

/**
 * @brief Function used to assign pointers to the total state density and 
 * the derivative
 * 
 * @param mu the chemical potential
 * @param p parameters for the equations - the temperature
 */
void NTotalAndDeriv(double mu, void* p, double *nTotal, double *nTotalDeriv)
{

    struct NTotalEvalParams * params = (struct NTotalEvalParams *)p;

    double t = (params->t);
    arma::mat energy = (params->energy);
    double nK = (params->nK); 
    double doping = (params->doping);

    *nTotal = 2.0*arma::accu(NFermi(mu, t, energy))/pow(nK, 2.0) - doping;
    *nTotalDeriv = 2.0*arma::accu(NFermiDeriv(mu, t, energy))/pow(nK, 2.0);
}

/**
 * @brief Function to append a value to a vector, note this is highly inefficient and should be
 * modified in the final code
 * 
 * @param vector vector to which the value is being appended to
 * @param value value being appended to the vector
 */
void vecPush(arma::vec& vector, const double& value) 
{

    arma::vec addValue(1);
    addValue.at(0) = value;
    vector.insert_rows(vector.n_rows, addValue.row(0));
    //addValue.print();

}

arma::cx_cube FftShift(arma::cx_cube& a, const int& dim)
{
    //get the length of the array in the dimension of the shift
    int length;

    //shifted output
    arma::cx_cube shiftA;
    shiftA.copy_size(a);

    //get real and imaginary parts
    arma::cube realA = arma::real(a);
    arma::cube imagA = arma::imag(a);
    
    if(dim == 1)
    {
        length = a.n_rows;

        //get the length of the shift
        int shiftMin = int(floor((double)length/2.0));

        //apply shift to each slice, this is done with lambda functions
        realA.each_slice([&shiftMin](arma::mat& tempA){tempA = shift(tempA, shiftMin, 0);});
        imagA.each_slice([&shiftMin](arma::mat& tempA){tempA = shift(tempA, shiftMin, 0);});

        arma::cx_cube shiftATemp(realA, imagA);
        
        shiftA = shiftATemp;
    }
    else if(dim == 2)
    {
        length = a.n_cols;

        //get the length of the shift
        int shiftMin = int(floor((double)length/2.0));

        //apply shift to each slice, this is done with lambda functions
        realA.each_slice([&shiftMin](arma::mat& tempA){tempA = shift(tempA, shiftMin, 1);});
        imagA.each_slice([&shiftMin](arma::mat& tempA){tempA = shift(tempA, shiftMin, 1);});

        arma::cx_cube shiftATemp(realA, imagA);
        
        shiftA = shiftATemp;
    }
    else if(dim == 3)
    {
        length = a.n_slices;

        //get the length of the shift
        int shiftMax = int(ceil((double)length/2.0));

        //get slices to swap
        arma::cx_cube tempA1 = a.slices(0, shiftMax);
        arma::cx_cube tempA2 = a.slices(shiftMax + 1, length - 1);

        //create new cube with swapped slices
        arma::cx_cube shiftATemp =  arma::join_slices(tempA1, tempA2);

        shiftA = shiftATemp;
    }

    return shiftA;

}

arma::cx_cube IfftShift(arma::cx_cube& a, const int& dim)
{
    //get the length of the array in the dimension of the shift
    int length;

    //shifted output
    arma::cx_cube shiftA;
    shiftA.copy_size(a);

    //get real and imaginary parts
    arma::cube realA = arma::real(a);
    arma::cube imagA = arma::imag(a);
    
    if(dim == 1)
    {
        length = a.n_rows;

        //get the length of the shift
        int shiftMin = int(floor((double)length/2.0));

        //apply shift to each slice, this is done with lambda functions
        realA.each_slice([&shiftMin](arma::mat& tempA){tempA = shift(tempA, -shiftMin, 0);});
        imagA.each_slice([&shiftMin](arma::mat& tempA){tempA = shift(tempA, -shiftMin, 0);});

        arma::cx_cube shiftATemp(realA, imagA);
        
        shiftA = shiftATemp;
    }
    else if(dim == 2)
    {
        length = a.n_cols;

        //get the length of the shift
        int shiftMin = int(floor((double)length/2.0));

        //apply shift to each slice, this is done with lambda functions
        realA.each_slice([&shiftMin](arma::mat& tempA){tempA = shift(tempA, -shiftMin, 1);});
        imagA.each_slice([&shiftMin](arma::mat& tempA){tempA = shift(tempA, -shiftMin, 1);});

        arma::cx_cube shiftATemp(realA, imagA);
        
        shiftA = shiftATemp;
    }
    else if(dim == 3)
    {
        length = a.n_slices;

        //get the length of the shift
        //get the length of the shift
        int shiftMin = int(floor((double)length/2.0));

        //get slices to swap
        arma::cx_cube tempA1 = a.slices(0, shiftMin);
        arma::cx_cube tempA2 = a.slices(shiftMin + 1, length - 1);

        //create new cube with swapped slices
        arma::cx_cube shiftATemp =  arma::join_slices(tempA1, tempA2);

        shiftA = shiftATemp;
    }

    return shiftA;
}

/**
 * @brief Function to pad a cube with slices of of value = val
 * 
 * @param a The cube being padded
 * @param nSlices The 'thickness' of the padding i.e. the number of slices being padded with
 * @param val The value being padded with
 * @return paddedA The padded cube
 */
arma::cx_cube PadCube(arma::cx_cube& a, const int& nSlices, const double& val)
{

    //create padding
    arma::cx_cube padding(a.n_rows, a.n_cols, nSlices);

    padding.fill(val);

    arma::cx_cube paddedA =  arma::join_slices(a, padding);

    return paddedA;


}


/**
 * @brief Construct a complex cube out of the real inptu data
 * 
 * @param in a real cube 
 * @return out a complex cube with imaginary part 0 and real part equal to in
 */
arma::cx_cube RealToComplex(const arma::cube& in)
{

    //create a cube of zeros of matching size
    arma::cube zeros;
    zeros.copy_size(in);

    //set all values to 0
    zeros.fill(0);

    arma::cx_cube out(in, zeros);
    
    return out;

}

/**
 * @brief Function to symmetrise a complex cube
 * 
 * @param in the complex cube to be symmetrised
 * @return out the symmetrised cube
 */
arma::cx_cube Symmetrise(arma::cx_cube& in)
{

    //create cube to be tranposed
    arma::cx_cube transposeIn;

    transposeIn = in;

    //tranpose each slice within the cube
    transposeIn.each_slice([](arma::cx_mat& tempA){tempA = tempA.st();});

    //calculate the symmetrised cube
    arma::cx_cube out = (in + transposeIn)/2.0;

    return out;

}


/**
 * @brief Function to symmetrise a cube
 * 
 * @param in the cube to be symmetrised
 * @return out the symmetrised cube
 */
arma::cube Symmetrise(arma::cube& in)
{

    //create cube to be tranposed
    arma::cube transposeIn;

    transposeIn = in;

    //tranpose each slice within the cube
    transposeIn.each_slice([](arma::mat& tempA){tempA = tempA.t();});

    //calculate the symmetrised cube
    arma::cube out = (in + transposeIn)/2;

    return out;

}

/**
 * @brief Function to transpose a cube
 * 
 * @param in the cube to be tranposed
 * @return tranpose the transposed cube
 */
arma::cube Transpose(arma::cube& in)
{

    //create cube to be tranposed
    arma::cube transpose = in;

    //tranpose each slice within the cube
    transpose.each_slice([](arma::mat& tempA){tempA = tempA.t();});

    return transpose;

}

/**
 * @brief Function to transpose a complex cube
 * 
 * @param in the cube to be tranposed
 * @return tranpose the transposed cube
 */
arma::cx_cube Transpose(arma::cx_cube& in)
{

    //create cube to be tranposed
    arma::cx_cube transpose = in;

    //tranpose each slice within the cube
    transpose.each_slice([](arma::cx_mat& tempA){tempA = tempA.st();});

    return transpose;

}

/**
 * @brief Function to clean diagonals of a cube
 * 
 * @param in input cube
 * @return out cleaned cube
 */
arma::cube CleanDiagonal(arma::cube& in)
{

    arma::cube out = in;

    //manually set diagonal elements to 0
    out.each_slice([](arma::mat& tempA)
    {

        arma::mat multip = arma::ones(size(tempA)) - arma::eye(size(tempA));
        tempA = tempA%multip;
        tempA = tempA%arma::flipud(multip);
    });

    return out;
}

/**
 * @brief Function to clean diagonals of a cube
 * 
 * @param in input cube
 * @return out cleaned cube
 */
arma::cx_cube CleanDiagonal(arma::cx_cube& in)
{

    arma::cx_cube out = in;

    //manually set diagonal elements to 0
    out.each_slice([](arma::cx_mat& tempA)
    {

        arma::mat multip = arma::ones(size(tempA)) - arma::eye(size(tempA));
        tempA = tempA%multip;
        tempA = tempA%arma::flipud(multip);
    });

    return out;
}

/**
 * @brief Flipud for a cube
 * 
 * @param in input cube
 * @return out flipud cube
 */
arma::cube FlipUDCube(arma::cube& in)
{

    arma::cube out = in;

    //flipud each slice
    out.each_slice([](arma::mat& tempA){tempA = arma::flipud(tempA);});

    return out;
}

/**
 * @brief Fliplr for a cube
 * 
 * @param in input cube
 * @return out fliplr cube
 */
arma::cube FlipLRCube(arma::cube& in)
{

    arma::cube out = in;

    //fliplr each slice
    out.each_slice([](arma::mat& tempA){tempA = arma::fliplr(tempA);});

    return out;
}

/**
 * @brief Flipud for a complex cube
 * 
 * @param in input cube
 * @return out flipud cube
 */
arma::cx_cube FlipUDCube(arma::cx_cube& in)
{

    arma::cx_cube out = in;

    //flipud each slice
    out.each_slice([](arma::cx_mat& tempA){tempA = arma::flipud(tempA);});

    return out;
}

/**
 * @brief Fliplr for a complex cube
 * 
 * @param in input cube
 * @return out fliplr cube
 */
arma::cx_cube FlipLRCube(arma::cx_cube& in)
{

    arma::cx_cube out = in;

    //fliplr each slice
    out.each_slice([](arma::cx_mat& tempA){tempA = arma::fliplr(tempA);});

    return out;
}

/**
 * @brief 3D linear interpolation for a cube object
 * 
 * @param x original x coordinates
 * @param y original y coordinates
 * @param z original z coordinates
 * @param in input gridded cube
 * @param xi x coordinates for interpolation
 * @param yi y coorindates for interpolation
 * @param zi z coordaintes for interpolation
 * @return interp the interpolated input matrix
 */
arma::cube Interpolate3D(const arma::vec& x, const arma::vec& y, const arma::vec& z, const arma::cube& in, const arma::vec& xi, const arma::vec& yi, const arma::vec zi)
{
    //set interpolated cube
    //temporary cube for in plane interpolation
    arma::cube interpTemp(xi.size(), yi.size(), z.size());
    //final interpolation
    arma::cube interp(xi.size(), yi.size(), zi.size());

    //only interpolate if necessary
    if(accu(abs(x - xi)) == 0 && accu(abs(y - yi)) == 0)
    {
        interpTemp = in;
    }
    else
    {
        //interpolate in the x-y plane
        for(unsigned int i = 0; i < z.size(); i++)
        {
            arma::interp2(x, y, in.slice(i), xi, yi, interpTemp.slice(i));
        }

    }

    //interpolate across z axis
    for(unsigned int i = 0; i < xi.size(); i++)
    {
        for(unsigned int j = 0; j < yi.size(); j++)
        {

            arma::vec initData = interpTemp.tube(i,j);

            //holder for interpolated data
            arma::vec interpData;

            //interpolate data
            arma::interp1(z, initData, zi, interpData);

            //linear extrapolation if the data is outside region
            /************************************
             * 
             * This needs to be made better so it actually catches all cases
             * 
             * ***********************************/
            if(interpData.has_nan() == true)
            {
                interpData[0] = 2.0*interpData[1] - interpData[2];
                interpData[interpData.size() - 1] = 2.0*interpData[interpData.size() - 2] - interpData[interpData.size() - 3];
            }

            interp.tube(i,j) = interpData;
        }
    }

    return interp;

}

/**
 * @brief A method to convert an armadillo vector to an alglib
 * 
 * @param vector A vector containing data to be transferrred
 * @return alglib::real_1d_array containing data in vector
 */
alglib::real_1d_array ConvertToAlglib(const arma::vec& vector)
{

    //get length of arrays
    unsigned int lengthOutput = vector.size();
    //set up output vector
    alglib::real_1d_array vectorOut;
    vectorOut.setlength(lengthOutput);

    //transfer data
    for(unsigned int i = 0; i < lengthOutput; i++)
    {
        vectorOut[i] = vector[i];
    }

    return vectorOut;

}

/**
 * @brief A method to convert to an armadillo vector from an 
 * alglib array
 * 
 * @param array Input alglib array
 * @return arma::vec  containing data in vector 
 */
arma::vec ConvertToArma(const alglib::real_1d_array& array)
{
    //get length of arrays
    unsigned int lengthOutput = array.length();
    //set up output vector
    arma::vec vectorOut;
    vectorOut.set_size(lengthOutput);

    //transfer data
    for(unsigned int i = 0; i < lengthOutput; i++)
    {
        vectorOut[i] = array[i];
    }

    return vectorOut;

}

/**
 * @brief Method to interpolate the lambda v temperature data to help find a better 
 * approximation for Tc. This interpolation function is a general interpolation function
 * which either performs linear interpolation through the inbuilt armadillo functions
 * or using the alglib cubic spline package.
 * 
 * @param vectorIn Input vector to be interpolated
 * @param vectorCoordsIn coordinates of the input vector
 * @param vectorOut output interpolated vector
 * @param vectorCoordsOut output coordinates of the interpolation
 * @param method linear or cubic
 */
void Interpolate1D(const arma::vec& vectorIn, const arma::vec& vectorCoordsIn, arma::vec& vectorOut, arma::vec& vectorCoordsOut, const std::string& method)
{

    if(method == "cubic")
    {
        //convert to alglib vectors for interpolation
        alglib::real_1d_array vectorInA = ConvertToAlglib(vectorIn);
        alglib::real_1d_array vectorCoordsInA = ConvertToAlglib(vectorCoordsIn);
        alglib::real_1d_array vectorOutA = ConvertToAlglib(vectorOut);
        alglib::real_1d_array vectorCoordsOutA = ConvertToAlglib(vectorCoordsOut);

        //calculate cubic spline
        alglib::spline1dconvcubic(vectorCoordsInA, vectorInA, vectorCoordsOutA, vectorOutA);

        //convert back to armadillo vector
        vectorOut = ConvertToArma(vectorOutA);
    }
    else if(method == "linear")
    {
        //interpolate the lambda data
        arma::interp1(vectorCoordsIn, vectorIn, vectorCoordsOut, vectorOut);

    }

}


/**
 * @brief Template to determine the sign of a number
 * 
 * @tparam T 
 * @param x
 * @return int 
 */
template <typename T> int sgn(T x) {
    return (T(0) < x) - (x < T(0));
}

#endif