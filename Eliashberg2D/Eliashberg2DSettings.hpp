/**
 * @file Eliashberg2DSettings.hpp
 * @author Theo Weinberger
 * @brief Code to read in settings for the Eliashberg code from a .cfg file
 * 
 * This file is passed a string from Eliashberg2DMain corresponding to a file 
 * name which contains the necessary data to run the Eliashberg script 
 * @version 0.1
 * @date 2021-11-09
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef ELIASHBERG2DSETTINGS_H
#define ELIASHBERG2DSETTINGS_H

#define ARMA_NO_DEBUG
#define ARMA_USE_MKL_ALLOC
#define ARMA_BLAS_LONG
#define ARMA_BLAS_LONG_LONG
#define ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_USE_MKL_TYPES
#define ARMA_ALLOW_FAKE_GCC

#include <iostream>
#include <mkl.h>
#include <armadillo>
#include <libconfig.h++>

/**
 * @brief Function to read in data that is to be stored in an armadillo vector
 * 
 * @param myVector Generic vector into which data should be put.
 * @param mySetting The configuration data that is to be read into myVector
 */
void ReadVector(arma::vec& myVector, const libconfig::Setting& mySetting);


/**
 * @brief Method to read file settings into Eliashberg code
 * 
 * @param fileName String data containing the name of the file containing configuration data.
 * @param t the tight binding hopping matrix element in meV 
 * @param ratioTight the ratio of tight-binding hopping elements tPrime/t
 * @param a The side length of a square lattice unit cell
 * @param doping The doping, calculated by Luttinger's therorem
 * @param magModel Whether the FM or AFM model is being used, takes values FM or AFM
 * @param phiModel The model being used for the anomalous self energy, takes values s, p, d
 * @param chi0 the static susceptibility
 * @param k0Squared k0^2 where k0 is the inverse correlation length without strong magnetic correlations
 * @param tSF The characteristic spin fluctuation temperature. Tsf*kappa02 is approximately constant
 * @param gSquaredChi0tSample g^2 chi0/t, can be input as an array
 * @param kSquaredSample k0^2 where k0 is the inverse correlation length without strong magnetic correlations
 * @param t0  Initial temperature for sampling
 * @param n0 The initial number of positive fermion Matsubara frequencies within the cutoff
 * @param omegaC The cutoff frequency, this should be around 20t but is equal to 2*pi*T0*N(T0)
 * @param nK The number of k points along kx, defining the amount of k-point sampling note at this is a 2D system the total number of points goes as Nk^2
 * @param errSigma The relative tolerance for the termination criteria for the convergence in calculating sigma
 * @param errLambda The relative tolerance for the termination criteria for the convergence in calculating lambda
 * @param maxIter The maximum number of iterations in finding self-consistent solutions at each temperature
 * @param phiRatio the ratio between scale of initial trial Phi and dPhi used in solving for lambda
 * @param symm whether sigma, the quasiparticle self energy, must be symmetrised
 * @param maxMuIter Max number of iterations when solving for mu
 * @param numLambdaSeg The number of segments that the final lambda data and corresponding temperature data are interpolated to
 * @param plot An array containing the mesh of the relaxation weight. This is scanned by the program in search of convergent solutions
 * @param relaxation  Whether we are plotting with respect to k or g2chi0t takes values of k or g
 * @return int Exit code
 */
int ReadFile(const std::string& fileName, double& t, double& ratioTight, double& a, double& doping, std::string& magModel, std::string& phiModel, double& chi0, double& k0Squared, double& tSF, arma::vec& gSquaredChi0tSample, arma::vec& kSquaredSample, double& t0, int& n0, double& omegaC, int& nK, double& errSigma, double& errLambda, int& maxIter, double& phiRatio, bool& symm, int& maxMuIter, int& numLambdaSeg, std::string& plot, arma::vec& relaxation);


#endif
