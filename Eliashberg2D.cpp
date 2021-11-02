/**
 * @file Eliashberg2D.cpp
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
#include "Eliashberg2D.hpp"


/**
 * @brief Construct a new Eliashberg object using the default paramters
 * 
 */
Eliashberg::Eliashberg()
{

    //Set default parameters for comparison against Ran's code
    //tight binding paramters
    _t = 1000/(4+8*0.45);
    _ratioTight = 0.45;
    _tPrime = _t*_ratioTight;
    _a = 1.0;

    //initialising parameters
    _t0 = 128*0.4*_t;
    _n0 = 8;
    _chi0 = 1.0;

    //cutoff frequency and sample
    _omegaC = 2*_n0*0.4*M_PI*_t;
     _nK = 16;

    //system physical parameters
    _tSF = 0.67*_t;
    _doping = 1.1;

    //sampling
    _gSquaredChi0tSample = 80.0;
    _k0Squared = 12.0;
    _kSquaredSample = 1.0; 

    //convergence control values
    _errSigma = 1e-10;
    _errLambda = 1e-15;
    _maxIter = 3000; 
    _relaxation = {0.1, 0.2, 0.5};
    _phiRatio = 1e8;
    _symm = true;
    _maxMuIter = 1000;

    //bools and vars controlling simulations
    _magModel = "FM";
    _phiModel = "d";
    _plot = "g";

}



void Eliashberg::SolveEliashberg()
{

    //initialise q space, 
    arma::vec qX = arma::linspace(-M_PI/_a*(1.0-1.0/_nK), M_PI/_a*(1.0-1.0/_nK), _nK);
    arma::vec qY = qX;

    _energy = _Dispersion(qX, qY);

    //initialise parameters for search
    double tC = 0.0;
    _kSquared = _kSquaredSample;
    double gSquaredChi0t = _gSquaredChi0tSample;
    double gSquared = gSquaredChi0t*_t/_chi0;

    double coupling; //coupling

    if(_magModel == "FM")
    {
        _qHat = _QMinus(qX, qY);
        coupling = gSquared/3;
    }
    else if(_magModel == "AFM")
    {
        _qHat = _QPlus(qX, qY);
        coupling = -gSquared;
    }

    /*initialise the relaxation weight index
    The relaxation weight smooths the transition between 
    steps to allow for more stable convergence to self-consistency*/
    unsigned int relaxIndex = 0;

    while(relaxIndex < _relaxation.size())
    {
        bool relaxIndexAccepted = true;
        /*the critical temperture Tc is defined where lambda = 1
        start the search for self consistent solutions where this holds*/
        double t = _t0;
        double n = _n0;

        //calculate the array containing the chemical potential
        arma::vec muArray = _calcMu(t, 0); //note 0 is picked as an initial trial value out of convenience

        /*Vectors containing the lambda values and temperature steps
        these start off as being unit length and are appended to each 
        iterations. Note the appending operation is highly inefficient 
        and so should be altered in later code*/
        arma::vec tStep;
        arma::vec lambdaVec(1, arma::fill::zeros);

        unsigned int muIndex = 0;

        //initialise the RG corrections 
        arma::cx_cube dSigma(_nK, _nK, 2*_n0, arma::fill::zeros); 
        arma::cube dPhi(_nK, _nK, 2*_n0, arma::fill::zeros); 

        /*Search for self consistency  of the Green's function
        is for when the value of lambda tends to 1. Therefore
        the search ends when lambda = 1*/
        while(abs(lambdaVec(lambdaVec.size() - 1)) < 1.0)
        {
            //intialise range of Matsubara frequencies
            arma::vec wMatsu = arma::linspace(-M_PI*(2*n - 1)*t, M_PI*(2*n - 1)*t, 2*n);
            arma::vec vMatsu = arma::linspace(-M_PI*(4*n - 2)*t, M_PI*(4*n - 2)*t, 4*n - 1);

            /* If the temperature, t, has droppped below the solvable region
            then use the most recent value for mu. This is because mu will converge
            as t tends to 0 and so the most recent value can be used as an approximate
            value for the real mu*/
            if(muIndex > muArray.size() - 1)
            {
                _mu = muArray(muArray.size() - 1);
            }
            else
            {
                _mu = muArray(muIndex);
            }

            //Intialise G
            arma::cx_cube gInv(_nK, _nK, _nK);
            arma::cx_cube g(_nK, _nK, _nK);

            for(int i = 0; i < _nK; i++)
            {
                for(int j = 0; j < _nK; j++)
                {
                    for(int k = 0; k < _nK; k++)
                    {
                        gInv(i,j,k) = (I*wMatsu(i) - (_energy(j,k) - _mu));
                    }
                }
            }

            g = 1/gInv;

            //Calculate chinu data, this corresponds to the analystical intergral of the dynamical suscpetibility
            arma::cube chiQNu =  _ChiQNu(vMatsu, qX, qY);

            //allocate memory for sigma
            arma::cube sigma(_nK, _nK, 2*_n0, arma::fill::zeros);

            //solve for sigma
            for(int i = 0; i < _maxIter; i++)
            {
                //evaluate new Greens function
                g = 1/(gInv - sigma);

                arma::cx_cube chiQNuComplex = RealToComplex(chiQNu);

                //Calculate convolution
                arma::cx_cube sigmaMatsuConv = gSquared*t/(pow(_nK, 2.0))*_MatsuConv(chiQNuComplex, g, 2*n, 4*n - 1) + dSigma;



            }

           



        }

    }












}


/**
 * @brief Function that calculates the energy of the system from its
 * dispersion relationsip
 * @param pX Momentum space vectors in the x direction
 * @param pY Momentum space vectors in the y direction
 * @return dispersion: a vector containing the dispersion relation of the system as a function of qX and qY
 */
arma::mat Eliashberg::_Dispersion(const arma::vec& pX, const arma::vec& pY)
{

    //solve dispersion relation for input q-space vectors, disersion relation is from tight-binding Hamiltonian

    //initialise the dispersion matrix
    arma::mat dispersion(_nK, _nK);

    for(int i = 0; i < _nK; i++)
    {

        for(int j = 0; j < _nK; j++)
        {
            dispersion(i, j) = -2.0*_t*(cos(pX[i]*_a) + cos(pY[j]*_a)) - 4.0*_tPrime*cos(pX[i]*_a)*cos(pY[j]*_a);
        }
    }


    return dispersion;

}

/**
 * @brief Function that calculates the negative root of momentum 
 * part of the susceptibility
 * 
 * @param qX Momentum space vectors in the x direction
 * @param qY Momentum space vectors in the y direction
 * @return qHatMinus: a vector containing the negative root of the susceptibility pole
 */
arma::vec Eliashberg::_QMinus(const arma::vec& qX, const arma::vec& qY)
{

    //calculate qHat, this applies in the FM case
    arma::vec qHatMinus = sqrt(4 - 2*(cos(qX*_a) + cos(qY*_a)));

    return qHatMinus;

}

/**
 * @brief Function that calculates the positive root of momentum 
 * part of the susceptibility
 * 
 * @param qX Momentum space vectors in the x direction
 * @param qY Momentum space vectors in the y direction
 * @return qHatPlus: a vector containing the positive root of the susceptibility pole
 */
arma::vec Eliashberg::_QPlus(const arma::vec& qX, const arma::vec& qY)
{

    //calculate qHat, this applies in the AFM case
    arma::vec qHatPlus = sqrt(4 + 2*(cos(qX*_a) + cos(qY*_a)));

    return qHatPlus;

}

/**
 * @brief Function that calculates the negative root of momentum 
 * part of the susceptibility. This version is for iterative calculations
 * 
 * @param qX Momentum space vectors in the x direction
 * @param qY Momentum space vectors in the y direction
 * @return qHatMinus: a vector containing the negative root of the susceptibility pole
 */
double Eliashberg::_QMinus(const double& qX, const double& qY)
{

    //calculate qHat, this applies in the FM case
    double qHatMinus = sqrt(4 - 2*(cos(qX*_a) + cos(qY*_a)));

    return qHatMinus;

}

/**
 * @brief Function that calculates the positive root of momentum 
 * part of the susceptibility. This version is for iterative calculations
 * 
 * @param qX Momentum space vectors in the x direction
 * @param qY Momentum space vectors in the y direction
 * @return qHatPlus: a vector containing the positive root of the susceptibility pole
 */
double Eliashberg::_QPlus(const double& qX, const double& qY)
{

    //calculate qHat, this applies in the AFM case
    double qHatPlus = sqrt(4 + 2*(cos(qX*_a) + cos(qY*_a)));

    return qHatPlus;

}

/**
 * @brief Function that calculates the eta parameter for the suscpetibility
 * eta is a momentum dependent coherence length (I think)
 * 
 * @param qX Momentum space vectors in the x direction
 * @param qY Momentum space vectors in the y direction
 */
void Eliashberg::_Eta(const arma::vec& qX, const arma::vec& qY)
{

    _eta = _tSF*_QMinus(qX, qY);

}

/**
 * @brief Function that calculates the eta parameter for the suscpetibility
 * eta is a momentum dependent coherence length (I think). This version 
 * is for iterative calculations
 * 
 * @param qX Momentum space vectors in the x direction
 * @param qY Momentum space vectors in the y direction
 */
double Eliashberg::_Eta(const double& qX, const double& qY)
{

    double eta = _tSF*_QMinus(qX, qY);

    return eta;

}

/**
 * @brief Function that calculates the cutoff frequency for the dynamics susceptibility integration
 * 
 * @return omega0: The cutoff frequency for the integration of the dynamic susceptibility
 */
arma::vec Eliashberg::_Omega0()
{
   
    arma::vec omega0 = _k0Squared*_eta;

    return omega0;

}

/**
 * @brief Function that calculates the cutoff frequency for the dynamics susceptibility integration
 * this version returns a double for iterative calucaltions
 * 
 * @param qX Momentum space vectors in the x direction
 * @param qY Momentum space vectors in the y direction
 * @return omega0: The cutoff frequency for the integration of the dynamic susceptibility
 */
double Eliashberg::_Omega0(const double& qX, const double& qY)
{
   
    double omega0 = _k0Squared*_Eta(qX, qY);

    return omega0;

}

/**
 * @brief Calculate the chemical potential energy by solving the eigenvalue equation
 * 
 * @param tInit The initial solver temperature
 * @param tTrial The trial solver temperature
 * 
 */
arma::vec Eliashberg::_calcMu(const double& tInit, const double& tTrial)
{
    //set up equation to solve
    double t = tInit;

    //ititialise storage vector for mu
    arma::vec muVec(1);

    //set status and counter
    int status = GSL_SUCCESS;
    int counter = 0;

    //root value
    double mu, muInit;

    //initialise roots
    mu = 1000;


    while(status == GSL_SUCCESS)
    {

        //define function to be solved
        gsl_function_fdf nTotal;
        struct NTotalEvalParams params = { t, _energy, _nK, _doping};

        nTotal.f = &NTotal;
        nTotal.df = &NTotalDeriv;
        nTotal.fdf = &NTotalAndDeriv;
        nTotal.params = &params;

        //define root solver, Here a Brent Solver is used
        const gsl_root_fdfsolver_type * T = gsl_root_fdfsolver_steffenson;
        gsl_root_fdfsolver * solver = gsl_root_fdfsolver_alloc (T);


        //initialise and run toor solver
        gsl_set_error_handler_off();
        gsl_root_fdfsolver_set(solver, &nTotal, mu);

        status = GSL_CONTINUE;

        for(int i = 0; i <= _maxMuIter && status == GSL_CONTINUE; ++i) 
        {
            /* iterate one step of the solver */
            status = gsl_root_fdfsolver_iterate(solver);
            if (status != GSL_SUCCESS)
            {
                break;
            }

            muInit = mu;

            /* get the solver's current best solution and bounds */
            mu = gsl_root_fdfsolver_root(solver);

            /* Check to see if the solution is within 0.001 */
            status = gsl_root_test_delta(mu, muInit, 0, 1e-10);

            if (status == GSL_SUCCESS)
            {

                //add mu to vector and reduce temperature

                //rather than appending to muVec, create a new vector at the desired length
                //and add values to this

                std::cout << mu << std::endl;
                arma::vec newMuVec(counter + 1);

                for(int j = 0; j < counter; j++)
                {
                    newMuVec[j] = muVec[j];
                }

                newMuVec[counter] = mu;

                muVec = newMuVec;                

                t = t/2.0;

            }

        }

        gsl_root_fdfsolver_free(solver);

        //index counter increment
        counter++;
    }

    return muVec;

}

/**
 * @brief Function to calculate the integral of the dynamic susceptibility
 * 
 * @param vmatsu The vector of matusbara frequencies
 * @param qx The mopmentum space vector in the x direction
 * @param qy The momentum space vector in the y direction
 * @return chiQNu The integrated dynamics susceptibility
 */
arma::cube Eliashberg::_ChiQNu(const arma::vec& vMatsu, const arma::vec& qX, const arma::vec& qY)
{
    arma::cube chiQNu(_nK, _nK, _nK);

    double omega0, eta, qHat, a;

    //Evaluate the integral of chiQNu at each point
    for(int i = 0; i < _nK; i++)
    {
        for(int j = 0; j < _nK; j++)
        {
            for(int k = 0; k < _nK; k++)
            {

                omega0 = _Omega0(qX[j], qY[k]);
                eta = _Eta(qX[j], qY[k]);

                //qhat value depends on magnetic model
                if(_magModel == "FM")
                {
                    qHat = _QMinus(qX[j], qY[k]);
                }
                else if(_magModel == "AFM")
                {
                    qHat = _QPlus(qX[j], qY[k]);
                }

                //last term in the integral of the dynamics susceptibility
                a = eta*(_kSquared + pow(qHat, 2.0));

                chiQNu(i, j, k) = (2*_chi0/M_PI)*_Omega0(qX[j], qY[k])*_ChiInt(omega0, vMatsu[i], a);

            }
        }
    }

    //replace NaN with 0
    chiQNu.replace(arma::datum::nan, 0);

    return chiQNu;

}

/**
 * @brief Function for the analytical expression of the dynamic susceptibility integral
 * 
 * @param y the variable that is being integrated over
 * @param a a system parameter
 * @param b a system parameter
 * @return The value of the integral at y
 */
double Eliashberg::_ChiInt(const double& y, const double& a, const double& b)
{

    double chiInt = (a*atan(y/a) - b*atan(y/b))/(pow(a, 2.0) - pow(b, 2.0));

    return chiInt;

}

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
 * @return matrixConv The matrix corresponding to a linear convolution in frequency of A and B
 */
arma::cube _MatsuConv(const arma::cx_cube& matrixA, const arma::cx_cube& matrixB, const int& lowindex, const int& highIndex)
{



    


}


int main()
{
    Eliashberg eliashberg;

    eliashberg.SolveEliashberg();

    return 0;
}



