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
    arma::mat _nFermi(const double& mu, const double& t);

    /**
     * @brief Function used to evaluate the total state density
     * 
     * @param mu the chemical potential
     * @param p parameters for the equations - the temperature
     * @return muSolve: the solved for value of mu
     */
    double _nTotalEval(double mu, void* );

    /**
     * @brief Function to solve an eigenvalue equation to find potential 
     * values of mu
     * 
     * @param tInit Initial temperature
     * @param doping The doping which corresponds to density of electron states
     * @param tTrial Current trial temperature for solving mu
     * @return muVec: A vector containing the values of mu
     */
    arma::vec _calcMu(const double& tInit, const double& tTrial);

    /**
     * @brief Function to calculate the integral of the dynamic susceptibility
     * 
     * @param vmatsu The vector of matusbara frequencies
     * @param qx The mopmentum space vector in the x direction
     * @param qy The momentum space vector in the y direction
     * @return chiQNu The integrated dynamics susceptibility
     */
    arma::cube _ChiQNu(const arma::vec& vmatsu, const arma::vec& qx, const arma::vec& qy);

    /**
     * @brief Function for the analytical expression of the dynamic susceptibility integral
     * 
     * @param y the variable that is being integrated over
     * @param a a system parameter
     * @param b a system parameter
     * @return The value of the integral at y
     */
    double _ChiInt(const double& y, const double& a, const double& b);

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
    double _gSquaredChi0tSample;

    /**
     * @brief k0Squared: k0^2 where k0 is the inverse correlation length without strong magnetic correlations
     * 
     */
    double _k0Squared;

    /**
     * @brief kSquaredSample: where k is the inverse correlation length with strong magnetic correlations, can be an array input
     * 
     */
    double _kSquaredSample;

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

arma::cx_cube FftShift(const arma::cx_cube& A, const int& dim)
{
    //get the length of the array in the dimension of the shift
    int length;

    //shifted output
    arma::cx_cube aShift;
    aShift.copy_size(A);
    
    if(dim == 1)
    {
        length = A.n_rows;

        //get the length of the shift
        int shiftMin = int(floor((double)length/2.0));
        int shiftMax = int(ceil((double)length/2.0));

        ////////////////////////////////////////////////////
        // iterate over slices and apply shift to each slice

        /////////////////////////////////////////////////////


    }
    else if(dim == 2)
    {
        length = A.n_cols;

        //get the length of the shift
        int shiftMin = int(floor((double)length/2.0));
        int shiftMax = int(ceil((double)length/2.0));
    }
    else if(dim == 3)
    {
        length = A.n_slices;

        //get the length of the shift
        int shiftMin = int(floor((double)length/2.0));
        int shiftMax = int(ceil((double)length/2.0));
    }




    return aShift;

}

arma::cx_cube IfftShift(const arma::cx_cube& A, const int& dim)
{
    //get the length of the array in the dimension of the shift
    int length;


    //shifted output
    arma::cx_cube aShift;
    aShift.copy_size(A);

    if(dim == 1)
    {
        length = A.n_rows;
    }
    else if(dim == 2)
    {
        length = A.n_cols;
    }
    else if(dim == 3)
    {
        length = A.n_slices;
    }

    //get the length of the shift
    int shiftLen = int(floor((double)length/2.0));


    return aShift;

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

#endif