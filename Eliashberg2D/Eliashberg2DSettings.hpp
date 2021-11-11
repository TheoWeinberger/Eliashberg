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

#define ARMA_NO_DEBUG

#include <iostream>
#include <armadillo>
#include <libconfig.h++>

#ifndef ELIASHBERG2DSETTINGS_H
#define ELIASHBERG2DSETTINGS_H

/**
 * @brief Function to read in data that is to be stored in an armadillo vector
 * 
 * @param myVector Generic vector into which data should be put.
 * @param mySetting The configuration data that is to be read into myVector
 */
void ReadVector(arma::vec& myVector, const libconfig::Setting& mySetting);

int ReadFile(const std::string& fileName, double& t, double& ratioTight, double& a, double& doping, std::string& magModel, std::string& phiModel, double& chi0, double& k0Squared, double& tSF, arma::vec& gSquaredChi0tSample, arma::vec& kSquaredSample, double& t0, int& n0, double& omegaC, int& nK, double& errSigma, double& errLambda, int& maxIter, double& phiRatio, bool& symm, int& maxMuIter, int& numLambdaSeg, std::string& plot, arma::vec& relaxation);


#endif
