/**
 * @file Eliashberg2DSettings.cpp
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

#include "Eliashberg2DSettings.hpp"

/**
 * @brief Function to read in data that is to be stored in an armadillo vector
 * 
 * @param myVector Generic vector into which data should be put.
 * @param mySetting The configuration data that is to be read into \a myVector
 */
void ReadVector(arma::vec& myVector, const libconfig::Setting& mySetting)
{
	int length = mySetting.getLength();
	
	myVector.set_size(length);

	for(int i = 0; i < length; i++)
	{
		double val  = double(mySetting[i]);
		myVector[i] = val;
	}
}


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
int ReadFile(const std::string& fileName, double& t, double& ratioTight, double& a, double& doping, std::string& magModel, std::string& phiModel, double& chi0, double& k0Squared, double& tSF, arma::vec& gSquaredChi0tSample, arma::vec& kSquaredSample, double& t0, int& n0, double& omegaC, int& nK, double& errSigma, double& errLambda, int& maxIter, double& phiRatio, bool& symm, int& maxMuIter, int& numLambdaSeg, std::string& plot, arma::vec& relaxation)
{
	//define storage variables
	libconfig::Config cfg;
	
	//read in diffraction configuration file, reporting error if it does not exist
	try
	{
	  	cfg.readFile(fileName.c_str());
	}
  	catch(const libconfig::FileIOException &fioex)
	{
		std::cerr << "I/O error while reading file." << std::endl;
		exit(EXIT_FAILURE);
	}
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
				<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}

  

	const libconfig::Setting& root = cfg.getRoot();

	//Get required data from file, these are values that must be set otherwise the code 
    //cannot run
	try
	{

        //tight binding parameters
        t = root["tightBinding"]["t"];
        ratioTight = root["tightBinding"]["ratioTight"];
        a = root["tightBinding"]["a"];
        doping = root["tightBinding"]["doping"];
        magModel = root["tightBinding"]["magModel"].c_str();
        phiModel = root["tightBinding"]["phiModel"].c_str();

        //general physical parameters
        chi0 = root["phyParam"]["chi0"];
        k0Squared = root["phyParam"]["k0Squared"];

        //get parameters being samples
        const libconfig::Setting& gSetting = root["sampleParam"].lookup("gSquaredChi0tSample");
		ReadVector(gSquaredChi0tSample, gSetting);

        const libconfig::Setting& kSetting = root["sampleParam"].lookup("kSquaredSample");
		ReadVector(kSquaredSample, kSetting);

        //get system sampling parameters
        n0 = root["samplingControl"]["n0"];
        nK = root["samplingControl"]["nK"];

        //get spin fluctuation temp
        //get convergence control parameters
        tSF = root["phyParam"]["tSF"];
        tSF *= t;

	}
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cerr << "Missing input value at: " << pex.getPath() << std::endl;
        exit(EXIT_FAILURE);
    } 

    //get convergence control parameters
    //convergence control parameters do not have to be set and can be left at default
    //get convergence control parameters
    //convergence control parameters do not have to be set and can be left at default
    try
    {
        //get convergence control parameters
        errSigma = root["convControl"]["errSigma"];
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
    } 

    //get next convergence control parameter
    try
    {
        //get convergence control parameters
        errLambda = root["convControl"]["errLambda"];
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
    } 

    //get next convergence control parameter
    try
    {
        //get convergence control parameters
        maxIter = root["convControl"]["maxIter"];
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
    } 

    //get next convergence control parameter
    try
    {
        //get convergence control parameters
        phiRatio = root["convControl"]["phiRatio"];
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
    } 

    //get next convergence control parameter
    try
    {
        //get convergence control parameters
        symm = root["convControl"]["symm"];
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
    } 

    //get next convergence control parameter
    try
    {
        //get convergence control parameters
        maxMuIter = root["convControl"]["maxMuIter"];
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
    } 

    //get next convergence control parameter
    try
    {
        //get convergence control parameters
        numLambdaSeg = root["convControl"]["numLambdaSeg"];
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
    }

    //get next convergence control parameter
    try
    {
        //get convergence control parameters
        const libconfig::Setting& relaxSetting = root["convControl"].lookup("relaxation");
		ReadVector(relaxation, relaxSetting);
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
    }

    //get remaining parameters, these can either be set in the cfg file 
    //however many of these are generally functions of other variables set earlier

    //get the initial temperature for sampling
    try
    {
        //get convergence control parameters
        t0 = root["samplingControl"]["t0"];
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
        t0 = 128*0.4*t;
    }

    //get the cutoff frequency for susceptibility integrations
    try
    {
        //get convergence control parameters
        omegaC = root["samplingControl"]["omegaC"];
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
        omegaC = 2*n0*0.4*M_PI*t;
    }

    //get plotting value
    try
    {
        //get convergence control parameters
        plot = root["output"]["g"].c_str();
    }
  	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
					<< " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
    catch(const libconfig::SettingNotFoundException &pex)
    {
        std::cout << "No input value at: " << pex.getPath() << std::endl;
        std::cout << "Setting to default" << std::endl;
    }

    return(EXIT_SUCCESS);
    
}
