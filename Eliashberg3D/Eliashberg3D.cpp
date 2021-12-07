/**
 * @file Eliashberg3D.cpp
 * @author Theo Weinberger
 * @brief This class uses the formalism from:
 * 
 *          Magnetically mediated superconductivity in quasi-two and three dimensions
 *          P. Monthoux and G. G. Lonzarich
 *          Phys. Rev. B 63, 054529 – Published 17 January 2001
 * 
 *          Magnetically mediated superconductivity: Crossover from cubic to tetragonal lattice
 *  *       P. Monthoux and G. G. Lonzarich
 *          Phys. Rev. B 66, 224504  – Published 12 December 2002
 * 
 * 
 * To calculate the energies, frequencies and transition temperatures in MeV and
 * the density of states in 1/meV for a spin mediated superconductor.
 * 
 * It is based of the code Tc_Adaptive.m developed by Ran Tao
 * 
 * This code requires a system to be reduced to the tight-binding approximation 
 * Hamiltonian with dispersion given by
 * 
 *   $\bm{\epsilon_p} = -2t[\cos(p_xa) + \cos(p_ya) + \alpha_t\cos(p_za)] - 4t'[\cos(p_xa)\cos(p_ya) 
 *    + \alpha_t\cos(p_xa)\cos(p_za) + \alpha_t\cos(p_ya)\cos(p_za)]$  
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
 * @date 2021-12-06
 * 
 * @copyright Copyright (c) 2021
 * 
 */

/*to do:

Find difference in dSigma

symmetrisationB

plotting/output


*/

#include "Eliashberg3D.hpp"
#include "Eliashberg3DSettings.hpp"


/**
 * @brief Construct a new Eliashberg object using the default paramters
 * 
 */
Eliashberg::Eliashberg()
{

    std::cout << "No configuration file given" << std::endl;
    std::cout << "Defaulting to example settings" << std::endl;

    //Set default parameters for comparison against Ran's code
    //tight binding paramters
    _t = 1000/(4+8*0.45);
    _ratioTight = 0.45;
    _tPrime = _t*_ratioTight;
    _a = 1.0;
    _alphaT = 1.0;

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
    _gSquaredChi0tSample = {80.0, 81.0, 82.0};
    _k0Squared = 12.0;
    _kSquaredSample = {1.0}; 

    //convergence control values
    _errSigma = 1e-5;
    _errLambda = 1e-5;
    _maxIter = 3000; 
    _relaxation = {0.1, 0.2, 0.5, 1.0};
    _phiRatio = 1e8;
    _symm = true;
    _maxMuIter = 1000;

    //How many steps to interpolate lambda and T to
    _numLambdaSeg = 1000000;

    //bools and vars controlling simulations
    _magModel = "AFM";
    _phiModel = "p";
    _plot = "g";

}


/**
 * @brief Construct a new Eliashberg object from input data from configuration 
 * file called, eliashberg.cfg
 * 
 */
Eliashberg::Eliashberg(const std::string& fileName)
{

    std::cout << "Beginning Eliashberg calculation from settings in " << fileName << std::endl;

    //Default convergence control values
    _errSigma = 1e-5;
    _errLambda = 1e-5;
    _maxIter = 3000; 
    _phiRatio = 1e8;
    _symm = true;
    _maxMuIter = 1000;
    _numLambdaSeg = 1000000;
    _relaxation = {0.1, 0.2, 0.5, 1.0};
    _alphaTSample = {1.0};
    _alphaMSample = {1.0};

    //set default plotting value
    _plot = "alpha";


    std::cout << " " << std::endl;

    ReadFile(fileName, _t, _ratioTight, _a, _alphaTSample, _alphaMSample, _doping, _magModel, _phiModel, _chi0,
    _k0Squared, _tSF, _gSquaredChi0tSample, _kSquaredSample, _t0, _n0, _omegaC, _nK,
    _errSigma, _errLambda, _maxIter, _phiRatio, _symm, _maxMuIter, _numLambdaSeg, _plot, _relaxation);

    //set tPrime according to ratioTight
    _tPrime = _t*_ratioTight;

    std::cout << " " << std::endl;


}


void Eliashberg::SolveEliashberg()
{

    //initialise q space, 
    arma::vec qX = arma::linspace(-M_PI/_a*(1.0-1.0/_nK), M_PI/_a*(1.0-1.0/_nK), _nK);
    arma::vec qY = qX;
    arma::vec qZ = qX;

    int lenKSample = _kSquaredSample.size();
    int lenGSample = _gSquaredChi0tSample.size();
    int lenAlphaT = _alphaTSample.size();
    int lenAlphaM = _alphaMSample.size();

    //check sampling is properly set up
    
    if(lenGSample > 1 && _plot == "k")
    {
        std::cout << "Data structure is not formatted for plotting in k space, please check the k sampling matrix and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of k, _plot should be set to k and the k sampling array should have multiple values" << std::endl;
        exit(1);
    }

    if(lenAlphaT > 1 && _plot == "k")
    {
        std::cout << "Data structure is not formatted for plotting in k space, please check the k sampling matrix and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of k, _plot should be set to k and the k sampling array should have multiple values" << std::endl;
        exit(1);
    }

    if(lenAlphaM > 1 && _plot == "k")
    {
        std::cout << "Data structure is not formatted for plotting in k space, please check the k sampling matrix and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of k, _plot should be set to k and the k sampling array should have multiple values" << std::endl;
        exit(1);
    }

    if(lenKSample > 1 && _plot == "g")
    {
        std::cout << "Data structure is not formatted for plotting in g space, please check the g sampling matrix and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of g, _plot should be set to g and the g sampling array should have multiple values" << std::endl;
        exit(1);
    }

    if(lenAlphaM > 1 && _plot == "g")
    {
        std::cout << "Data structure is not formatted for plotting in g space, please check the g sampling matrix and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of g, _plot should be set to g and the g sampling array should have multiple values" << std::endl;
        exit(1);
    }
    
    if(lenAlphaT > 1 && _plot == "g")
    {
        std::cout << "Data structure is not formatted for plotting in g space, please check the g sampling matrix and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of g, _plot should be set to g and the g sampling array should have multiple values" << std::endl;
        exit(1);
    }

    if(_plot == "kg" && lenAlphaT > 1)
    {
        std::cout << "Data structure is not formatted for plotting in k-g space, please check the sampling matrices and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of k-g, _plot should be set to kg and the alpha sampling arrays should be single valued" << std::endl;
        exit(1);
    }

    if(_plot == "kg" && lenAlphaM > 1)
    {
        std::cout << "Data structure is not formatted for plotting in k-g space, please check the sampling matrices and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of k-g, _plot should be set to kg and the alpha sampling arrays should be single valued" << std::endl;
        exit(1);       
    }

    if(_plot == "alpha" && lenGSample > 1)
    {
        std::cout << "Data structure is not formatted for plotting in alpha space, please check the sampling matrices and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of alpha, _plot should be set to alpha and the g and k sampling arrays should be single valued" << std::endl;
        exit(1);
    }

    if(_plot == "alpha" && lenKSample > 1)
    {
        std::cout << "Data structure is not formatted for plotting in alpha space, please check the sampling matrices and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of alpha, _plot should be set to alpha and the g and k sampling arrays should be single valued" << std::endl;
        exit(1);       
    }

    arma::mat tC, tInit;
    //temperature arrays
    if(_plot == "alpha")
    {
        tC.resize(lenAlphaM, lenAlphaT);
        tC.fill(0.0);
        /*********************************/
        tInit.resize(lenAlphaM, lenAlphaT);
        tInit.fill(0.0); //does this get used
    }
    else
    {
        tC.resize(lenKSample, lenGSample);
        tC.fill(0.0);
        /*********************************/
        tInit.resize(lenKSample, lenGSample);
        tInit.fill(0.0); //does this get used
    }

    //initialise energy
    _energy = _Dispersion(qX, qY, qZ);
    /**********************************/

    //calculate the array containing the chemical potential
    arma::vec muArray = _calcMu(_t0); //note 0 is picked as an initial trial value out of convenience

    //initialise remainder of variables now for more efficient run
    //physical/sampling parameters
    double gSquaredChi0t, gSquared, coupling;
    //temperature and intial number of matsubara frequencies
    double t, n;

    //counting parameters 
    unsigned int relaxIndex; //counting through the relaxation indeces
    int counter; //counter for the number of lambda iterations
    unsigned int muIndex; //counter for the chemical potential

    //allocate memory for sigma and phi
    std::vector<arma::cx_cube> sigma(2*_n0);
    std::vector<arma::cube> phi(2*_n0);
    std::vector<arma::cube> phiFilter(2*_n0);

    //intiialise matsubara convolutions
    std::vector<arma::cx_cube> sigmaMatsuConv(2*_n0);
    std::vector<arma::cube> phiMatsuConv(2*_n0);

    //factors in matsubara convolutions
    std::vector<arma::cx_cube> gAbsSquaredComp(2*_n0);
    std::vector<arma::cube> gAbsSquared(2*_n0);
    std::vector<arma::cx_cube> gAbsPhi(2*_n0);

    //initialise the RG corrections 
    std::vector<arma::cx_cube> dSigma(2*_n0); 
    std::vector<arma::cube> dPhi(2*_n0); 
    std::vector<arma::cube> dPhiTemp(2*_n0); 
    std::vector<arma::cube> dSigmaReal(2*_n0);
    std::vector<arma::cube> dSigmaImag(2*_n0);
    std::vector<arma::cube> dSigmaRealTemp(2*_n0);
    std::vector<arma::cube> dSigmaImagTemp(2*_n0);


    //Intialise G
    std::vector<arma::cx_cube> gInv(2*_n0);
    std::vector<arma::cx_cube> g(2*_n0);

    //temporary holder for MatsuConv 
    std::vector<arma::cx_cube> matsuConvTemp;


    //set size of cubes in each instance
    for(int i = 0; i < 2*_n0; i++)
    {
        sigma[i].set_size(_nK, _nK, _nK);
        phi[i].set_size(_nK, _nK, _nK);
        phiFilter[i].set_size(_nK, _nK, _nK);
        sigmaMatsuConv[i].set_size(_nK, _nK, _nK);
        phiMatsuConv[i].set_size(_nK, _nK, _nK);
        gAbsSquaredComp[i].set_size(_nK, _nK, _nK);
        gAbsSquared[i].set_size(_nK, _nK, _nK);
        gAbsPhi[i].set_size(_nK, _nK, _nK);
        dSigma[i].set_size(_nK, _nK, _nK);
        dPhi[i].set_size(_nK, _nK, _nK);
        dPhiTemp[i].set_size(_nK, _nK, _nK);
        dSigmaReal[i].set_size(_nK, _nK, _nK);
        dSigmaImag[i].set_size(_nK, _nK, _nK);
        dSigmaRealTemp[i].set_size(_nK, _nK, _nK);
        dSigmaImagTemp[i].set_size(_nK, _nK, _nK);
        gInv[i].set_size(_nK, _nK, _nK);
        g[i].set_size(_nK, _nK, _nK);
    }

    
    /*Vectors containing the lambda values and temperature steps
    these start off as being unit length and are appended to each 
    iterations. Note the appending operation is highly inefficient 
    and so should be altered in later code*/
    arma::vec tStep;
    arma::vec lambdaVec;

    //direct test against Ran code
    muArray = {1.852935461551102e+02,108.578641736149,94.6635609987612,110.255073698169,131.715764552851,145.752034210909,151.043635672090,152.067024877566,152.404983516244,152.193764044987,151.640222096965};

    for(int a = 0; a < lenKSample; a++)
    {
        for(int b = 0; b < lenGSample; b++)
        {
            for(int c = 0; c < lenAlphaT; c++)
            {
                for(int d = 0; d < lenAlphaM; d++)
                {
                    //set assymetry
                    _alphaT = _alphaTSample[c];
                    _alphaM = _alphaTSample[d];

                    if(_alphaT == 1.0 && _alphaM == 1.0)
                    {
                        _symmType = 2;
                    }
                    else
                    {
                        _symmType = 1;
                    }
              
                    //user information
                    std::cout << "Initialising run with parameters" << std::endl;
                    std::cout << "g: " << _gSquaredChi0tSample[b] << std::endl;
                    std::cout << "k: " << _kSquaredSample[a] << std::endl;

                    //initialise energy
                    _energy = _Dispersion(qX, qY, qZ);


                    //initialise parameters for search
                    _kSquared = _kSquaredSample[a];
                    gSquaredChi0t = _gSquaredChi0tSample[b];
                    gSquared = gSquaredChi0t*_t/_chi0;

                    if(_magModel == "FM")
                    {
                        coupling = gSquared/3;
                    }
                    else if(_magModel == "AFM")
                    {
                        coupling = -gSquared;
                    }
                    else
                    {
                        std::cout << "Invalid value of the magnetic model" << std::endl;
                        std::cout << "Please enter either FM or AFM" << std::endl;
                        exit(1);
                    }

                    /*initialise the relaxation weight index
                    The relaxation weight smooths the transition between 
                    steps to allow for more stable convergence to self-consistency*/
                    relaxIndex = 0;


                    while(relaxIndex < _relaxation.size())
                    {

                        std::cout << " "  << std::endl;
                        std::cout << "Attempt with relaxation value of: " << _relaxation[relaxIndex] << std::endl; 
                        std::cout << " "  << std::endl;

                        //set size of tStep and lambdaVec to 1
                        tStep.set_size(1);
                        lambdaVec.set_size(1);

                        bool relaxIndexAccepted = true;
                        /*the critical temperture Tc is defined where lambda = 1
                        start the search for self consistent solutions where this holds*/

                        //set the inital temperature
                        t = _t0;
                        //and set initial number of matsubara frequencies
                        n = _n0;

                        //set the counter for the chemical potential index to zero
                        muIndex = 0;

                        //initialise the RG corrections 
                        for(int i = 0; i < 2*_n0; i++)
                        {
                            dSigma[i].zeros();
                            dPhi[i].zeros();
                        }

                        //counter for iteration number
                        counter = 0;

                        //boolean to see if it is the first run
                        bool firstRun;

                        //boolean to state whether regular matsubara convolution is occuring or the RG corrections
                        bool rg = false;

                        /*Search for self consistency  of the Green's function
                        is for when the value of lambda tends to 1. Therefore
                        the search ends when lambda = 1. The loop is broken
                        if counter goes above 40 where Tc ~ 0*/
                        while(abs(lambdaVec(lambdaVec.size() - 1)) < 1.0  && counter < 40)
                        {

                            //reset all cubes to be filled with 0
                            for(int i = 0; i < 2*_n0; i++)
                            {
                                sigma[i].zeros();
                                phi[i].zeros();
                                phiFilter[i].zeros();
                            }

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

                            //note this order of looping is optimal
                            for(int l = 0; l < 2*_n0; l++)
                            {
                                for(int k = 0; k < _nK; k++)
                                {
                                    for(int j = 0; j < _nK; j++)
                                    {
                                        for(int i = 0; i < _nK; i++)
                                        {
                                            //define complex double 
                                            std::complex<double> val(- (_energy(i,j,k) - _mu), wMatsu(l));
                                            gInv[l](i,j,k) = val;
                                        }
                                    }
                                }
                            }

                            for(int l = 0; l< 2*_n0; l++)
                            {
                                g[l] = 1.0/gInv[l];
                            }

                            //Calculate chinu data, this corresponds to the analystical intergral of the dynamical suscpetibility
                            std::vector<arma::cube> chiQNu =  _ChiQNu(vMatsu, qX, qY, qZ);
                            std::vector<arma::cx_cube> chiQNuComplex = RealToComplex(chiQNu);

                            
                            //double to store the relative error in the convergence in sigma
                            double relErrS = INFINITY;

                            //solve for sigma
                            for(int i = 0; i < _maxIter; i++)
                            {
                                //evaluate new Greens function
                                for(int l = 0; l < 2*_n0; l++)
                                {
                                    g[l] = 1/(gInv[l] - sigma[l]);
                                }

                                //set firsrRun to true if it is the first run. This signals that the fftw plans should be set
                                if(counter == 0 && i == 0)
                                {
                                    firstRun = true;
                                }

                                //Calculate convolution
                                /************************************************************************************************
                                 *  Weird thing with incompatible arrays here, check with Ran
                                 * There may also be an error in the greens function
                                 */
                                matsuConvTemp = _MatsuConv(chiQNuComplex, g, 2*n, 4*n - 1, firstRun, rg);
                                for(int l = 0; l < 2*_n0; l++)
                                {
                                    sigmaMatsuConv[l] = gSquared*t/(pow(_nK, 3.0))*matsuConvTemp[l] + dSigma[l];
                                }

                                //set firstRun to false for remaining iterations
                                firstRun = false;

                                //work out relative erorr in convergence
                                double deltaS = 0;
                                double totalS = 0;
                                
                                for(int l = 0; l < 2*_n0; l++)
                                {
                                    deltaS += arma::accu(abs(sigmaMatsuConv[l]-sigma[l]));
                                    totalS += arma::accu(abs(sigma[l]));
                                }
                                relErrS = deltaS/totalS;

                                //iterate the sigma matrix
                                //switched relaxation because this makes more sense
                                for(int l = 0; l < 2*_n0; l++)
                                {
                                    sigma[l] = (1 - _relaxation[relaxIndex])*sigmaMatsuConv[l] + _relaxation[relaxIndex]*sigma[l]; 
                                }                    
                            
                                /*sigma should be symmetric in rows and columns (i.e. within the slices)
                                make sure that this is the case as it can vary a bit over time due to 
                                floating point/rounding errors*/
                                if(_symm == true)
                                {
                                    if(_symmType == 1)
                                    {
                                        sigma = SymmetriseA(sigma);
                                    }
                                    else if(_symmType == 2)
                                    {
                                        sigma = SymmetriseB(sigma);
                                    }
                                    else
                                    {
                                        std::cout << "Unknown input for symmetrisation type of sigma matrix" << std::endl;
                                        std::cout << "Please choose either symmetrisation type A or B" << std::endl;
                                        exit(1);
                                    }
                                }          

                                if(relErrS < _errSigma)
                                {
                                    break;
                                }

                            }
                            
                            /**************************************************************/
                            //this does not seem to make sense and I think it should send the loop back to the beginning already
                            //i.e.. a break should be inserted in the else statement
                            if(relErrS > _errSigma)
                            {

                                //exit failure if no convergence can be achieved
                                if (relaxIndex == _relaxation.size() - 1)
                                {
                                    std::cout << "No convergent relaxation weight could be found within the array of relaxation values. Please provide a larger maxIter, a finer mesh of relxation weights, or a larger relative error tolerance to achieve convergence" << std::endl;
                                    exit(1);
                                }
                                else
                                {

                                    relaxIndexAccepted = false;
                                    relaxIndex += 1;
                                    break;
                                }

                            }


                            //find the `eigenvalue' lambda using the power methid
                            //set the phi depending on the phi model
                            phi = _PhiFun(qX, qY);
                            phiFilter = _PhiSymm(qX, qY);

                

                            //clean diagonals
                            if(_phiModel == "d")
                            {

                                phi = CleanDiagonal(phi);
                                phiFilter =  CleanDiagonal(phiFilter);

                            }

                            //scale phi
                            arma::vec maxTemp(2*_n0);

                            for( int l = 0; l < 2*_n0; l++)
                            {
                                maxTemp[l] = abs(dPhi[l]).max();
                            }
                            double scale = abs(maxTemp).max();

                            //scale phi if needed
                            if(scale > 0)
                            {
                                for(int l = 0; l < 2*_n0; l++)
                                {
                                    phi[l] *= scale*_phiRatio;
                                }
                            }

                            //calcualte <Phi|Phi> , shoudln;t this be a proper inner product???
                            /*************************************************************/
                            double phiInner = 0;
                            for(int l = 0; l < 2*_n0; l++)
                            {
                                phiInner += arma::accu(phi[l]%phi[l]);
                            }

                            //get the absolute value of g
                            for(int l = 0; l < 2*_n0; l++)
                            {
                                gAbsSquared[l] = pow(abs(g[l]), 2.0);
                            }

                            gAbsSquaredComp = RealToComplex(gAbsSquared);

                            //convolve chi wth |G|^2%phi
                            for(int l = 0; l < 2*_n0; l++)
                            {
                                gAbsPhi[l] = gAbsSquaredComp[l]%phi[l];
                            }

                            /********************
                             * 
                             * in general phimatsu conv can be imaginary, change this
                             **********************************/
                            matsuConvTemp = _MatsuConv(chiQNuComplex, gAbsPhi, 2*n, 4*n - 1, firstRun, rg);
                            for(int l = 0; l < 2*_n0; l++)
                            {
                                phiMatsuConv[l] = arma::real(coupling*t/pow(_nK, 3.0)*matsuConvTemp[l]) + dPhi[l];
                            }

                            double lambda = 0;
                            //calcilate lambda as <Phi|Phi> = <Phi|A|Phi> = <Phi|Phi1>
                            for(int l = 0; l < 2*_n0; l++)
                            {
                                lambda += arma::accu(phi[l]%phiMatsuConv[l])/phiInner;
                            }

                            //apply symmetry transforms
                            phi = _SymmByFiltLabel(phiMatsuConv, phiFilter);

                            //initialise relative error and normalisation for lambda search
                            double relErrL = INFINITY;
                            bool normalise = false;

                            for(int i = 0; i < _maxIter; i++)
                            {

                                phiInner = 0;

                                for(int l = 0; l < 2*_n0; l++)
                                {
                                    //calcualte <Phi|Phi>
                                    phiInner += arma::accu(phi[l]%phi[l]);
                                }

                                for(int l = 0; l < 2*_n0; l++)
                                {
                                    //convolve chi wth |G|^2%phi
                                    gAbsPhi[l] = gAbsSquaredComp[l]%phi[l];
                                }

                                matsuConvTemp = _MatsuConv(chiQNuComplex, gAbsPhi, 2*n, 4*n - 1, firstRun, rg);

                                for(int l = 0; l < 2*_n0; l++)
                                {
                                    phiMatsuConv[l] = arma::real(coupling*t/pow(_nK, 3.0)*matsuConvTemp[l]) + dPhi[l];
                                }

                                //calcilate lambda as <Phi|Phi> = <Phi|A|Phi> = <Phi|Phi1>
                                double lambda1 = 0;

                                for(int l = 0; l < 2*_n0; l++)
                                {
                                    lambda1 += arma::accu(abs(phi[l]%phiMatsuConv[l]))/phiInner;
                                }

                                //calculate the relative change in lambda
                                double relErrL1 = abs(lambda1 - lambda)/abs(lambda);
                                
                                //check normalisation cases to update lambda
                                if(normalise == true)
                                {
                                    relErrL = relErrL1;

                                    //iterate to the next step
                                    phi = _SymmByFiltLabel(phiMatsuConv, phiFilter);
                                    for(int l = 0; l < 2*_n0; l++)
                                    {
                                        phi[l] /= lambda1;
                                    }

                                    lambda = lambda1;
                                }
                                else if(relErrL1 > relErrL)
                                {
                                    //if relative error increases, set to normalise
                                    normalise = true;
                                    for(int l = 0; l < 2*_n0; l++)
                                    {
                                        phi[l] = phi[l]/lambda;
                                    }
                                }
                                else
                                {
                                    //unormalised iteration
                                    relErrL = relErrL1;
                                    phi = _SymmByFiltLabel(phiMatsuConv, phiFilter);
                                    lambda = lambda1;
                                }

                                //break loop if threshold reached
                                if(relErrL < _errLambda)
                                {
                                    break;
                                }

                            }

                            if(relErrL > _errLambda)
                            {
                                std::cout << "No convergent wavefunction phi could be found. Please provide a larger maxIter, a larger ratio, or a larger relative error tolerance to achieve convergence" << std::endl;
                                exit(1);
                            }
                            //normalise if not done alread
                            else if(normalise == false)
                            {
                                for(int l = 0; l < 2*_n0; l++)
                                {
                                    phi[l] /= lambda;
                                }
                            }
                            
                            //rather than appending to lambdaVec, create a new vector at the desired length
                            //and add values to this
                            arma::vec newLambdaVec(counter + 1);

                            for(int j = 0; j < counter; j++)
                            {
                                newLambdaVec[j] = lambdaVec[j];
                            }

                            std::cout << "For iteration " << counter << " lambda: " << lambda << std::endl;

                            newLambdaVec[counter] = lambda;

                            lambdaVec = newLambdaVec;  

                            /* Now calculate the RG corrections. Since T is halved at the 
                            next step we also know the mesh of fermionic frequencies is 
                            halved too*/
                            int m = int(ceil(n/2.0));

                            //get the fermionic frequencies in the L domain

                            arma::vec wL = wMatsu.rows((n - m),(n + m - 1));

                            //cut mastubara v frequencies
                            arma::vec vCut = arma::linspace(-(4*m - 2)*M_PI*t, (4*m - 2)*M_PI*t, 4*m - 1);

                            //Calculate chinu data, this corresponds to the analystical intergral of the dynamical suscpetibility
                            //this is now calculated in the cut frequency domain
                            std::vector<arma::cube> chiCut =  _ChiQNu(vCut, qX, qY, qZ);
                            std::vector<arma::cx_cube> chiCutComplex = RealToComplex(chiCut);

                            //crop g 
                            std::vector<arma::cx_cube> gL(g.begin() + (n - m), g.begin() + (n + m));

                            rg = true;

                            //calculate the contribution to this iteration from L to L domains
                            matsuConvTemp = _MatsuConv(chiCutComplex, gL, 2*m, 4*m - 1, firstRun, rg);

                            std::vector<arma::cx_cube> sigmaLL(matsuConvTemp.size());
                            for(unsigned int l = 0; l < sigmaLL.size(); l++)
                            {
                                sigmaLL[l] = gSquared*t/(pow(_nK, 3.0))*matsuConvTemp[l];
                            }

                            //contribution to sigma from L to H domains
                            std::vector<arma::cx_cube> sigmaL(sigma.begin() + (n - m), sigma.begin() + (n + m)); 


                            //holder for dSigma

                            //note I think 3D interpolation is unecessary as it is just along the frequency domain that interpolation occurs
                            //interpolate real and imaginary parts
                            for(unsigned int l = 0; l < sigmaL.size(); l++)
                            {
                                arma::cx_cube diffTemp = sigmaL[l] - sigmaLL[l];
                                dSigmaRealTemp[l] = arma::real(diffTemp);
                                dSigmaImagTemp[l] = arma::imag(diffTemp);
                            }

                            dSigmaReal = Interpolate4D(qX, qY, qZ, wL, dSigmaRealTemp, qX, qY, qZ, wMatsu/2.0, "cubic");
                            dSigmaImag = Interpolate4D(qX, qY, qZ, wL, dSigmaImagTemp, qX, qY, qZ, wMatsu/2.0, "cubic");

                            //set whole dSigma matrix
                            for(unsigned int l = 0; l < dSigmaImag.size(); l++)
                            {
                                dSigma[l].set_real(dSigmaReal[l]);
                                dSigma[l].set_imag(dSigmaImag[l]);
                            }

                            //perform same cutting for phi

                            std::vector<arma::cube> phiL(phi.begin() + (n - m), phi.begin() + (n + m)); 

                            std::vector<arma::cx_cube> gAbsPhiL(phiL.size());
                            
                            for(unsigned int l = 0; l < phiL.size(); l++)
                            {
                                gAbsPhiL[l] = RealToComplex(pow(abs(gL[l]), 2.0))%phiL[l];
                            }

                            matsuConvTemp = _MatsuConv(chiCutComplex, gAbsPhiL, 2*m, 4*m - 1, firstRun, rg);

                            std::vector<arma::cube> phiLL(matsuConvTemp.size());
                            for(unsigned int l = 0; l < matsuConvTemp.size(); l++)
                            {
                                phiLL[l] = arma::real(coupling*t/pow(_nK, 3.0)*matsuConvTemp[l]);
                            }

                            for(unsigned int l = 0; l < phiLL.size(); l++)
                            {
                                dPhiTemp[l] = lambda*phiL[l] - phiLL[l];
                            }

                            dPhi = Interpolate4D(qX, qY, qZ, wL, dPhiTemp, qX, qY, qZ, wMatsu/2.0, "cubic");

                            //rather than appending to tStep, create a new vector at the desired length
                            //and add values to this
                            arma::vec newTStep(counter + 1);

                            for(int j = 0; j < counter; j++)
                            {
                                newTStep[j] = tStep[j];
                            }

                            newTStep[counter] = t;

                            tStep = newTStep;

                            t /= 2.0;

                            //step the mu index
                            muIndex++;

                            counter++;
                        }
                        
                            
                        //break the loop if this relaxation index is suitable
                        if(relaxIndexAccepted == true)
                        {
                            break;
                        }

                    }

                    //query points for the temperature
                    arma::vec tInterp = arma::linspace(2*t, 4*t, 1000000);
                    arma::vec lInterp;
                    lInterp.copy_size(tInterp);
                    lInterp.zeros();

                    //interpolate the lambda data
                    Interpolate1D(lambdaVec, tStep, lInterp, tInterp, "cubic");

                    //Tc the index where lQuery first goes above 1
                    arma::uword indexTC = (abs(lInterp - 1.0)).index_min();

                    //extract the critical temperature
                    if(_plot == "alpha")
                    {
                        tC(d, c) = tInterp(indexTC);
                        //output the critical temperature
                        std::cout << "    " << std::endl;
                        std::cout << "Tc: "<< tC(d,c) << std::endl;
                        std::cout << "Tc/Tsf: " << tC(d,c)/_tSF << std::endl;
                        std::cout << "    " << std::endl;
                    }
                    else
                    {
                        tC(a, b) = tInterp(indexTC);
                        //output the critical temperature
                        std::cout << "    " << std::endl;
                        std::cout << "Tc: "<< tC(a,b) << std::endl;
                        std::cout << "Tc/Tsf: " << tC(a,b)/_tSF << std::endl;
                        std::cout << "    " << std::endl;
                    }
                }
            }
        }      
    } 

    //output k/g data if relevant
    if(lenKSample > 1 && _plot == "k")
    {
        std::cout << "k: " << std::endl;
        _kSquaredSample.t().print();

        //output all the critical temperature
        std::cout << "Tc: " << std::endl;
        tC.t().print();
        std::cout << "Tc/Tsf: " << std::endl;
        (tC/_tSF).t().print();
    }

    if(lenGSample > 1 && _plot == "g")
    {
        std::cout << "g: " << std::endl;
        _gSquaredChi0tSample.t().print();

        //output all the critical temperature
        std::cout << "Tc: " << std::endl;
        tC.print();
        std::cout << "Tc/Tsf: " << std::endl;
        (tC/_tSF).print();
    }



    //output data to csv file for plotting
    //plotting as a function of gSquared
    if(_plot == "g")
    {
        arma::mat outputData(lenGSample, 2);   

        for(int i = 0; i < lenGSample; i++)
        {

            outputData(i, 0) = _gSquaredChi0tSample[i];
            outputData(i, 1) = tC[i]/_tSF;
            
        }

        outputData.save("gData_k" + std::to_string(_kSquaredSample[0]) + "_" + _magModel, arma::csv_ascii);
    }
    //plot as a function of kappa
    else if(_plot == "k")
    {
        arma::mat outputData(lenKSample, 2);   

        for(int i = 0; i < lenKSample; i++)
        {

            outputData(i, 0) = _kSquaredSample[i];
            outputData(i, 1) = tC[i]/_tSF;
            
        }

        outputData.save("kData_g" + std::to_string(_gSquaredChi0tSample[0]) + "_" + _magModel, arma::csv_ascii);
    }
    else if(_plot == "kg")
    {
        arma::mat outputData = _ScaleTC(tC, _kSquaredSample, _gSquaredChi0tSample);

        outputData.save("kgData_" + _magModel, arma::csv_ascii);
    }
    else if(_plot == "alpha")
    {
        arma::mat outputData = _ScaleTC(tC, _alphaMSample, _alphaTSample);

        outputData.save("alphaData_" + _magModel, arma::csv_ascii);

    }
    else
    {
        std::cout << "Incorrect plotting values, please set plot to either g, k, kg or alpha" << std::endl;
        exit(1);
    }

    _DeleteDFTPlans();

}


/**
 * @brief Function that calculates the energy of the system from its
 * dispersion relationsip
 * @param pX Momentum space vectors in the x direction
 * @param pY Momentum space vectors in the y direction
 * @param pZ Momentum space vectors in the z direction
 * @return dispersion: a vector containing the dispersion relation of the system as a function of qX and qY
 */
arma::cube Eliashberg::_Dispersion(const arma::vec& pX, const arma::vec& pY, const arma::vec& pZ)
{

    //solve dispersion relation for input q-space vectors, disersion relation is from tight-binding Hamiltonian

    //initialise the dispersion matrix
    arma::cube dispersion(_nK, _nK, _nK);

    for(int k = 0; k < _nK; k++)
    {    
        for(int j = 0; j < _nK; j++)
        {
            for(int i = 0; i < _nK; i++)
            {
                dispersion(i, j, k) = -2.0*_t*(cos(pX[i]*_a) + cos(pY[j]*_a) + _alphaT*cos(pZ[k]*_a)) - 4.0*_tPrime*(cos(pX[i]*_a)*cos(pY[j]*_a) + _alphaT*cos(pX[i]*_a)*cos(pZ[k]*_a) + _alphaT*cos(pZ[k]*_a)*cos(pY[j]*_a));
            }
        }
    }


    return dispersion;

}


/**
 * @brief Function that calculates the negative root of momentum 
 * part of the susceptibility. This version is for iterative calculations
 * 
 * @param qX Momentum space vectors in the x direction
 * @param qY Momentum space vectors in the y direction
 * @return qHatMinus: a vector containing the negative root of the susceptibility pole
 */
double Eliashberg::_QMinus(const double& qX, const double& qY, const double& qZ)
{

    //calculate qHat, this applies in the FM case
    double qHatMinus = sqrt((4 + 2*_alphaM) - 2*(cos(qX*_a) + cos(qY*_a) + _alphaM*cos(qZ*_a)));

    return qHatMinus;

}


/**
 * @brief Function that calculates the positive root of momentum 
 * part of the susceptibility. This version is for iterative calculations
 * 
 * @param qX Momentum space vector in the x direction
 * @param qY Momentum space vector in the y direction
 * @param qZ Momentum space vector in the z direction
 * @return qHatPlus: a vector containing the positive root of the susceptibility pole
 */
double Eliashberg::_QPlus(const double& qX, const double& qY, const double& qZ)
{

    //calculate qHat, this applies in the AFM case
    double qHatPlus = sqrt((4 + 2*_alphaM) + 2*(cos(qX*_a) + cos(qY*_a) + _alphaM*cos(qZ*_a)));

    return qHatPlus;

}


/**
 * @brief Function that calculates the eta parameter for the suscpetibility
 * eta is a momentum dependent coherence length (I think)/. 
 * This version is for iterative calculations
 * 
 * @param qX Momentum space vector in the x direction
 * @param qY Momentum space vector in the y direction
 * @param qZ Momentum space vector in the z direction
 * @return eta: a vector containing the coherence length values
 */
double Eliashberg::_Eta(const double& qX, const double& qY, const double& qZ)
{

    double eta = _tSF*_QMinus(qX, qY, qZ);

    return eta;

}


/**
 * @brief Function that calculates the cutoff frequency for the dynamics susceptibility integration.
 * This version is for iterative calculations
 * 
 * @param qX Momentum space vector in the x direction
 * @param qY Momentum space vector in the y direction
 * @param qZ Momentum space vector in the z direction
 * @return omega0: The cutoff frequency for the integration of the dynamic susceptibility
 */
double Eliashberg::_Omega0(const double& qX, const double& qY, const double& qZ)
{
   
    double omega0 = _k0Squared*_Eta(qX, qY, qZ);

    return omega0;

}


/**
 * @brief Calculate the chemical potential energy by solving the eigenvalue equation
 * 
 * @param tInit The initial solver temperature
 */
arma::vec Eliashberg::_calcMu(const double& tInit)
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
    mu = 100;


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
            status = gsl_root_test_delta(mu, muInit, 0, 1e-1);

            if (status == GSL_SUCCESS)
            {

                //add mu to vector and reduce temperature

                //rather than appending to muVec, create a new vector at the desired length
                //and add values to this
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
 * @param qZ Momentum space vector in the z direction
 * @return chiQNu The integrated dynamics susceptibility
 */
std::vector<arma::cube> Eliashberg::_ChiQNu(const arma::vec& vMatsu, const arma::vec& qX, const arma::vec& qY, const arma::vec& qZ)
{

    int iMax = qX.n_elem;
    int jMax = qY.n_elem;
    int kMax = qZ.n_elem;
    int lMax = vMatsu.n_elem;

    std::vector<arma::cube> chiQNu(lMax);

    for(int i = 0; i < lMax; i++)
    {
        chiQNu[i].set_size(iMax, jMax, kMax);
    }

    double omega0, eta, qHat, a;
    
    //Evaluate the integral of chiQNu at each point
    for(int l = 0; l < lMax; l++)
    {
        for(int k = 0; k < kMax; k++)
        {
            for(int j = 0; j < jMax; j++)
            {
                for(int i = 0; i < iMax; i++)
                {

                    omega0 = _Omega0(qX[i], qY[j], qZ[k]);
                    eta = _Eta(qX[i], qY[j], qZ[k]);

                    //qhat value depends on magnetic model
                    if(_magModel == "FM")
                    {
                        qHat = _QMinus(qX[i], qY[j], qZ[k]);
                    }
                    else if(_magModel == "AFM")
                    {
                        qHat = _QPlus(qX[i], qY[j], qZ[k]);
                    }
                    else
                    {
                        std::cout << "Invalid value of the magnetic model" << std::endl;
                        std::cout << "Please enter either FM or AFM" << std::endl;
                        exit(1);
                    }

                    //last term in the integral of the dynamics susceptibility
                    a = eta*(_kSquared + pow(qHat, 2.0));

                    chiQNu[l](i, j, k) = (2*_chi0/M_PI)*omega0*_ChiInt(omega0, vMatsu[l], a);

                }
            }
        }
    }

    //replace NaN with 0
    for(int i = 0; i < lMax; i++)
    {
        chiQNu[i].replace(arma::datum::nan, 0);
    }

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
 * @param firstRun boolean specifying if it is the first
 * @param boolean specifying if rg corrections are occuring 
 * @return matrixConv The matrix corresponding to a linear convolution in frequency of A and B
 */
std::vector<arma::cx_cube> Eliashberg::_MatsuConv(const std::vector<arma::cx_cube>& matrixA, const std::vector<arma::cx_cube>& matrixB, const int& lowIndex, const int& highIndex, const bool& firstRun, const bool& rg)
{

    //These tempoarary cubes will store the padded data 
    std::vector<arma::cx_cube> matrixAPadded = matrixA;
    std::vector<arma::cx_cube> matrixBPadded = matrixB;

    //Apply circular shifts to get frequency centred output

    for(unsigned int i = 0; i < matrixA.size(); i ++)
    {
        matrixAPadded[i] = IfftShift(matrixAPadded[i], 1);
        matrixAPadded[i] = IfftShift(matrixAPadded[i], 2);
        matrixAPadded[i] = IfftShift(matrixAPadded[i], 3);
    }

    for(unsigned int i = 0; i < matrixB.size(); i ++)
    {
        matrixBPadded[i] = IfftShift(matrixBPadded[i], 1);
        matrixBPadded[i] = IfftShift(matrixBPadded[i], 2); 
        matrixBPadded[i] = IfftShift(matrixBPadded[i], 3); 
    }

    //pad matrices
    matrixAPadded = Pad(matrixAPadded, matrixB.size() - 1, 0.0);
    matrixBPadded = Pad(matrixBPadded, matrixA.size() - 1, 0.0);

    arma::cx_vec convolution;

    if(rg == true)
    {

        int nElemPlan = matrixAPadded.size()*matrixAPadded[0].n_elem;
        //cube for the planner to use and overwrite
        arma::cx_vec planner(nElemPlan);

        int n[] {int(matrixAPadded.size()), int(matrixAPadded[0].n_slices), int(matrixAPadded[0].n_cols), int(matrixAPadded[0].n_rows)};

        _SetDFTPlansRG(planner, planner, n);

        //flatten 
        arma::cx_vec matrixAFlat = Flatten(matrixAPadded);
        arma::cx_vec matrixBFlat = Flatten(matrixBPadded);

        //apply ffts, multiple and then apply ifft to achieve convolution
        fftw_execute_dft(_forwardPlanRG, (double(*)[2])&matrixAFlat(0), (double(*)[2])&matrixAFlat(0));
        fftw_execute_dft(_forwardPlanRG, (double(*)[2])&matrixBFlat(0), (double(*)[2])&matrixBFlat(0));

        convolution = matrixAFlat%matrixBFlat;

        fftw_execute_dft(_inversePlanRG, (double(*)[2])&convolution(0), (double(*)[2])&convolution(0));
    }
    else
    {

        //set DFTs if it is the first run
        if(firstRun == true)
        {
            int nElemPlan = matrixAPadded.size()*matrixAPadded[0].n_elem;
            //cube for the planner to use and overwrite
            arma::cx_vec planner(nElemPlan);

            int n[] {int(matrixAPadded.size()), int(matrixAPadded[0].n_slices), int(matrixAPadded[0].n_cols), int(matrixAPadded[0].n_rows)};

            _SetDFTPlans(planner, planner, n);
        }

        //flatten 
        arma::cx_vec matrixAFlat = Flatten(matrixAPadded);
        arma::cx_vec matrixBFlat = Flatten(matrixBPadded);

        //apply ffts, multiple and then apply ifft to achieve convolution
        fftw_execute_dft(_forwardPlan, (double(*)[2])&matrixAFlat(0), (double(*)[2])&matrixAFlat(0));
        fftw_execute_dft(_forwardPlan, (double(*)[2])&matrixBFlat(0), (double(*)[2])&matrixBFlat(0));

        convolution = matrixAFlat%matrixBFlat;

        fftw_execute_dft(_inversePlan, (double(*)[2])&convolution(0), (double(*)[2])&convolution(0));

    }

    convolution /= convolution.n_elem;

    //truncate the frequency domain to remain within the cutoff range
    std::vector<arma::cx_cube> convolutionReshape = Make4D(convolution, matrixAPadded);
    
    std::vector<arma::cx_cube> convolutionCrop(convolutionReshape.begin() + lowIndex - 1, convolutionReshape.begin() + highIndex);

    //reshift matrix
    for(unsigned int i = 0; i < convolutionCrop.size(); i++)
    {
        convolutionCrop[i] = IfftShift(convolutionCrop[i], 1);
        convolutionCrop[i] = IfftShift(convolutionCrop[i], 2);
        convolutionCrop[i] = IfftShift(convolutionCrop[i], 3);
    }    

    return convolutionCrop; 
}


/**
 * @brief Function to set DFT plans for the matsubara frequency convolution
 * 
 * @param in The matrix being transformed
 * @param out The output matrix of the DFT
 * @param n array containing dimesions
 */
void Eliashberg::_SetDFTPlans(const arma::cx_vec& in, const arma::cx_vec& out, const int n[])
{

    _forwardPlan = fftw_plan_dft(4, n, (double(*)[2])&in(0), (double(*)[2])&out(0), FFTW_FORWARD, FFTW_MEASURE);

    _inversePlan = fftw_plan_dft(4, n, (double(*)[2])&in(0), (double(*)[2])&out(0), FFTW_BACKWARD, FFTW_MEASURE);

}


/**
 * @brief Function to set DFT plans for the matsubara frequency convolution
 * 
 * @param in The matrix being transformed
 * @param out The output matrix of the DFT
 * @param n array containing dimesions
 */
void Eliashberg::_SetDFTPlansRG(const arma::cx_vec& in, const arma::cx_vec& out, const int n[])
{


    _forwardPlanRG = fftw_plan_dft(4, n, (double(*)[2])&in(0), (double(*)[2])&out(0), FFTW_FORWARD, FFTW_ESTIMATE);

    _inversePlanRG = fftw_plan_dft(4, n, (double(*)[2])&in(0), (double(*)[2])&out(0), FFTW_BACKWARD, FFTW_ESTIMATE);

}


/**
 * @brief Function to delete DFT plans for the matsubara frequency convolution
 * 
 */
void Eliashberg::_DeleteDFTPlans()
{

    fftw_destroy_plan(_forwardPlan);

    fftw_destroy_plan(_inversePlan);

    fftw_destroy_plan(_forwardPlanRG);

    fftw_destroy_plan(_inversePlanRG);

}


/**
 * @brief Calculate the phi function used to iteratively calculate the eigenalue lambda
 * via the power method
 * 
 * @param qX A vector containing the momentum in the x direction
 * @param qY A vector conaining the momentum in the y direction
 * @return phi A cube containing the data for phi
 */
std::vector<arma::cube> Eliashberg::_PhiFun(const arma::vec& qX, const arma::vec& qY)
{

    std::vector<arma::cube> phi(2*_n0);

    for(int i = 0; i < 2*_n0; i++)
    {
        phi[i].set_size(_nK, _nK, _nK);
    }

    //case for s wave
    if(_phiModel == "s")
    {
        for(int i = 0; i < 2*_n0; i++)
        {
            phi[i].fill(1.0);
        }
    }
    //case for p wave
    else if(_phiModel == "p")
    {
        for(int l = 0; l < 2*_n0; l++)
        {
            for(int i = 0; i < _nK; i++)
            {
                for(int j = 0; j < _nK; j++)
                {
                    for(int k = 0; k < _nK; k++)
                    {

                        phi[l](k, j, i) = sin(qX[k]*_a);
                        
                    }
                }
            }    
        }    
    }
    //case for d wave
    else if(_phiModel == "d")
    {
        for(int l = 0; l < 2*_n0; l++)
        {
            for(int i = 0; i < _nK; i++)
            {
                for(int j = 0; j < _nK; j++)
                {
                    for(int k = 0; k < _nK; k++)
                    {

                        phi[l](k, j, i) = cos(qX[k]*_a) -  cos(qY[j]*_a);
                        
                    }
                }
            }   
        }
    }
    //catch errors
    else
    {
        std::cout << "Incorrect input for the orbital model, phiModel." << std::endl;
        std::cout << "Please choose a value of either s, p, or d." << std::endl;
        exit(1);
    }

    return phi;

}


/**
 * @brief Calculate the symmetric phi function used to iteratively calculate the eigenalue lambda
 * via the power method
 * 
 * @param qX A vector containing the momentum in the x direction
 * @param qY A vector conaining the momentum in the y direction
 * @return phi A cube containing the data for symmetric phi
 */
std::vector<arma::cube> Eliashberg::_PhiSymm(const arma::vec& qX, const arma::vec& qY)
{

    std::vector<arma::cube> phi(2*_n0);

    for(int i = 0; i < 2*_n0; i++)
    {
        phi[i].set_size(_nK, _nK, _nK);
    }

    //case for s wave
    if(_phiModel == "s")
    {
        for(int i = 0; i < 2*_n0; i++)
        {
            phi[i].fill(1.0);
        }
    }
    //case for p wave
    else if(_phiModel == "p")
    {
        for(int l = 0; l < 2*_n0; l++)
        {
            for(int i = 0; i < _nK; i++)
            {
                for(int j = 0; j < _nK; j++)
                {
                    for(int k = 0; k < _nK; k++)
                    {

                        phi[l](k, j, i) = sgn(qX[k]);
                        
                    }
                }
            } 
        }       
    }
    //case for d wave
    else if(_phiModel == "d")
    {
        for(int l = 0; l < 2*_n0; l++)
        {
            for(int i = 0; i < 2*_n0; i++)
            {
                for(int j = 0; j < _nK; j++)
                {
                    for(int k = 0; k < _nK; k++)
                    {

                        phi[l](k, j, i) = sgn(qX[k] - qY[j])*sgn(qX[k] + qY[j]);
                        
                    }
                }
            }  
        } 
    }
    //catch errors
    else
    {
        std::cout << "Incorrect input for the orbital model, phiModel." << std::endl;
        std::cout << "Please choose a value of either s, p, or d." << std::endl;
        exit(1);
    }

    return phi;

}


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
std::vector<arma::cube> Eliashberg::_SymmByFiltLabel(std::vector<arma::cube>& matrixA, const std::vector<arma::cube>& filter)
{

    //initialise symmetrised matrix
    int length = matrixA.size();
    std::vector<arma::cube> matrixS(length);

    //apply symmetry transforms depending on phimodel
    if(_phiModel == "s")
    {
        matrixS = matrixA;
    }
    else if(_phiModel == "p")
    {
        for(int i = 0; i < length; i++)
        {
            arma::cube matrixA1 = matrixA[i]%filter[i];
            arma::cube matrixA2 = FlipUDCube(matrixA1);
            matrixS[i] = (matrixA1 + matrixA2)%filter[i]/2.0;
        }
    }
    else if(_phiModel == "d")
    {
        for(int i = 0; i < length; i++)
        {
            arma::cube matrixA1 = matrixA[i]%filter[i];
            arma::cube matrixA2 = FlipLRCube(matrixA1);
            arma::cube matrixA3 = FlipUDCube(matrixA1);
            arma::cube matrixA4 = FlipLRCube(matrixA3);
            arma::cube matrixA5 = Transpose(matrixA1);
            arma::cube matrixA6 = FlipLRCube(matrixA5);
            arma::cube matrixA7 = FlipUDCube(matrixA5);
            arma::cube matrixA8 = FlipLRCube(matrixA7);
            matrixS[i] = (matrixA1 + matrixA2 + matrixA3 + matrixA4 + matrixA5 + matrixA6 + matrixA7 + matrixA8)%filter[i]/8.0;
        }
    }
    else
    {
        std::cout << "Incorrect input for the orbital model, phiModel." << std::endl;
        std::cout << "Please choose a value of either s, p, or d." << std::endl;
        exit(1);
    }

    return matrixS;

}

/**
 * @brief Add scales to temperature output
 * 
 * @param in input temperature
 * @param xScale scale for the rows
 * @param yScale scale for the columns
 * @return scaledMat: the temperature matrix with scales for plotting
 */
arma::mat Eliashberg::_ScaleTC(const arma::mat& in, const arma::vec& xScale, const arma::vec& yScale)
{
    arma::mat scaledMat = in;

    //normalise
    scaledMat /= _tSF;
    
    //generate containers to store axis data, width is offset by one to keep matrix rectangular
    int width = xScale.size();

    arma::rowvec xScaleAxis;
    xScaleAxis.set_size(width + 1);


    //shift axes that they are centred on zero
    for(int i = 0; i < width + 1; i++)
    {
        if(i == 0)
        {
            xScaleAxis[i] = 0;
        }
        else
        {
            xScaleAxis[i] = xScale[i - 1];
        }
    }

    //append to grating matrix. This is slow, must look into faster methods. 
    scaledMat.insert_cols(0, yScale);
    scaledMat.insert_rows(0, xScaleAxis);

    return scaledMat;
}


/**
 * @brief Function to solve for the number of fermions in the system
 * 
 * @param mu the chemical potential
 * @param t the temperature
 * @return arma::mat nFermiMat matrix of the nFermi at each sampling point
 */
arma::cube NFermi(const double& mu, const double& t, const arma::cube& energy)
{

    arma::cube nFermiMat = 1.0/(exp((energy-mu)/t) + 1);

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
arma::cube NFermiDeriv(const double& mu, const double& t, const arma::cube& energy)
{

    arma::cube nFermiMatDeriv = exp((energy-mu)/t)/(t*pow((exp((energy-mu)/t) + 1), 3.0));

    return nFermiMatDeriv;

}


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
    arma::cube energy = (params->energy);
    double nK = (params->nK); 
    double doping = (params->doping);

    double nTotal = 2.0*arma::accu(NFermi(mu, t, energy))/pow(nK, 3.0) - doping;

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
    arma::cube energy = (params->energy);
    double nK = (params->nK); 

    double nTotalDeriv = 2.0*arma::accu(NFermiDeriv(mu, t, energy))/pow(nK, 3.0);

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
    arma::cube energy = (params->energy);
    double nK = (params->nK); 
    double doping = (params->doping);

    *nTotal = 2.0*arma::accu(NFermi(mu, t, energy))/pow(nK, 3.0) - doping;
    *nTotalDeriv = 2.0*arma::accu(NFermiDeriv(mu, t, energy))/pow(nK, 3.0);
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

}


/**
 * @brief Apply fftshift in the dimension of the cube 
 * as specified. This is to centre the 0 frequency output 
 * of a fourier transform. The fft shift is equivalent to
 * a circular permutation of floor(a.length/2).
 * 
 * @param a The input cube
 * @param dim Dimension which the shift is being applied, 1 - rows, 2 - columns, 3 - slices.
 * @return arma::cx_cube Shifted cube
 */
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
        arma::cx_cube tempA1 = a.slices(0, shiftMax - 1);
        arma::cx_cube tempA2 = a.slices(shiftMax, length - 1);

        //create new cube with swapped slices
        arma::cx_cube shiftATemp =  arma::join_slices(tempA2, tempA1);

        shiftA = shiftATemp;
    }

    return shiftA;

}


/**
 * @brief Apply Ifftshift in the dimension of the cube 
 * as specified. This is to centre the 0 frequency output 
 * of a fourier transform. The fft shift is equivalent to
 * a circular permutation of - floor(a.length/2). Ifftshift will
 * always undo the operation of Fftshift
 * 
 * @param a The input cube
 * @param dim Dimension which the shift is being applied, 1 - rows, 2 - columns, 3 - slices.
 * @return arma::cx_cube Shifted cube
 */
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
        arma::cx_cube tempA1 = a.slices(0, shiftMin - 1);
        arma::cx_cube tempA2 = a.slices(shiftMin, length - 1);

        //create new cube with swapped slices
        arma::cx_cube shiftATemp =  arma::join_slices(tempA2, tempA1);

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
 * @brief Function to pad a cube with slices of of value = val
 * 
 * @param a The cube being padded
 * @param n The 'thickness' of the padding i.e. the number cubes being padded with
 * @param val The value being padded with
 * @return paddedA The padded cube
 */
std::vector<arma::cx_cube> Pad(std::vector<arma::cx_cube>& a, const int& n, const double& val)
{

    //create padding

    arma::cx_cube padding(a[0].n_rows, a[0].n_cols, a[0].n_slices);
    padding.fill(val);

    std::vector<arma::cx_cube> paddingVec(n, padding);

    std::vector<arma::cx_cube> paddedA;

    paddedA.reserve(a.size() + paddingVec.size()); //preallocate memory

    //concatenate
    paddedA.insert(paddedA.end(), a.begin(), a.end());
    paddedA.insert(paddedA.end(), paddingVec.begin(), paddingVec.end());

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
 * @brief Construct a complex cube out of the real inptu data
 * 
 * @param in a real cube 
 * @return out a complex cube with imaginary part 0 and real part equal to in
 */
std::vector<arma::cx_cube> RealToComplex(const std::vector<arma::cube>& in)
{

    //create a cube of zeros of matching size
    int length = in.size();
    std::vector<arma::cx_cube> out(length);

    arma::cube zeros;
    zeros.copy_size(in[0]);

    //set all values to 0
    zeros.fill(0);

    for(int i = 0; i < length; i++)
    {
        arma::cx_cube outTemp(in[i], zeros);
        
        out[i] = outTemp;
    }
    
    return out;

}


/**
 * @brief Function to symmetrise a complex cube
 * 
 * @param in the complex cube to be symmetrised
 * @return out the symmetrised cube
 */
std::vector<arma::cx_cube> SymmetriseA(std::vector<arma::cx_cube>& in)
{

    //create cube to be tranposed
    int length = in.size();

    std::vector<arma::cx_cube> transposeIn;
    std::vector<arma::cx_cube> out(length);

    transposeIn = in;

    for(int i = 0; i < length; i++)
    {
        //tranpose each slice within the cube
        transposeIn[i].each_slice([](arma::cx_mat& tempA){tempA = tempA.st();});

        //calculate the symmetrised cube
        out[i] = (in[i] + transposeIn[i])/2.0;

    }

    return out;

}


/**
 * @brief Function to symmetrise a cube
 * 
 * @param in the cube to be symmetrised
 * @return out the symmetrised cube
 */
std::vector<arma::cube> SymmetriseA(std::vector<arma::cube>& in)
{

    //create cube to be tranposed
    int length = in.size();

    std::vector<arma::cube> transposeIn;
    std::vector<arma::cube> out(length);

    transposeIn = in;

    for(int i = 0; i < length; i++)
    {
        //tranpose each slice within the cube
        transposeIn[i].each_slice([](arma::mat& tempA){tempA = tempA.t();});

        //calculate the symmetrised cube
        out[i] = (in[i] + transposeIn[i])/2.0;
    }

    return out;

}

/**
 * @brief Function to symmetrise a cube
 * 
 * @param in the cube to be symmetrised
 * @return out the symmetrised cube
 */
std::vector<arma::cx_cube> SymmetriseB(std::vector<arma::cx_cube>& in)
{
    //create cube to be tranposed
    int length = in.size();

    arma::cx_cube symm1;
    arma::cx_cube symm2;
    arma::cx_cube symm3;
    arma::cx_cube symm4;
    arma::cx_cube symm5;
    std::vector<arma::cx_cube> out(length);

    for(int i = 0; i < length; i++)
    {
        symm1 = Permute(in[i], "213");
        symm2 = Permute(in[i], "231");
        symm3 = Permute(in[i], "312");
        symm4 = Permute(in[i], "321");
        symm5 = Permute(in[i], "132");
        out[i] = (in[i] + symm1 + symm2 + symm3 + symm4 + symm5)/6.0; 
    }

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
std::vector<arma::cube> CleanDiagonal(std::vector<arma::cube>& in)
{

    std::vector<arma::cube> out = in;

    for(unsigned int i = 0; i < out.size(); i++)
    {
        //manually set diagonal elements to 0
        out[i].each_slice([](arma::mat& tempA)
        {

            arma::mat multip = arma::ones(size(tempA)) - arma::eye(size(tempA));
            tempA = tempA%multip;
            tempA = tempA%arma::flipud(multip);
        });
    }

    return out;
}


/**
 * @brief Function to clean diagonals of a cube
 * 
 * @param in input cube
 * @return out cleaned cube
 */
std::vector<arma::cx_cube> CleanDiagonal(std::vector<arma::cx_cube>& in)
{

    std::vector<arma::cx_cube> out = in;

    for(unsigned int i = 0; i < out.size(); i++)
    {
        //manually set diagonal elements to 0
        out[i].each_slice([](arma::cx_mat& tempA)
        {

            arma::mat multip = arma::ones(size(tempA)) - arma::eye(size(tempA));
            tempA = tempA%multip;
            tempA = tempA%arma::flipud(multip);
        });
    }

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
void Interpolate1D(const arma::vec& vectorIn, const arma::vec& vectorCoordsIn, arma::vec& vectorOut, const arma::vec& vectorCoordsOut, const std::string& method)
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
 * @brief 3D linear interpolation for a cube object
 * 
 * @param x original x coordinates
 * @param y original y coordinates
 * @param z original z coordinates
 * @param in input gridded cube
 * @param xi x coordinates for interpolation
 * @param yi y coorindates for interpolation
 * @param zi z coordaintes for interpolation
 * @param method cubic or linear
 * @return interp the interpolated input matrix
 */
arma::cube Interpolate3D(const arma::vec& x, const arma::vec& y, const arma::vec& z, const arma::cube& in, const arma::vec& xi, const arma::vec& yi, const arma::vec& zi, const std::string& method)
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
    for(unsigned int j = 0; j < yi.size(); j++)
    {
        for(unsigned int i = 0; i < yi.size(); i++)
        {

            arma::vec initData = interpTemp.tube(i,j);

            //holder for interpolated data
            arma::vec interpData;

            //interpolate data
            //arma::interp1(z, initData, zi, interpData);
            Interpolate1D(initData, z, interpData, zi, method);

            //linear extrapolation if the data is outside region
            if(interpData.has_nan() == true)
            {
                arma::uvec nonFiniteCoords = arma::find_nonfinite(interpData);

                if(nonFiniteCoords[0] == 0 && nonFiniteCoords[nonFiniteCoords.size() - 1] == interpData.size() - 1)
                {
                    interpData[0] = 2.0*interpData[1] - interpData[2];
                    interpData[interpData.size() - 1] = 2.0*interpData[interpData.size() - 2] - interpData[interpData.size() - 3];
                }
            }
            else if(interpData[0] == 0 && interpData[interpData.size() - 1] == 0)
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
 * @brief 4D linear interpolation for a cube object
 * 
 * @param x original x coordinates
 * @param y original y coordinates
 * @param z original z coordinates
 * @param w original frequency coordinates
 * @param in input gridded cube
 * @param xi x coordinates for interpolation
 * @param yi y coorindates for interpolation
 * @param zi z coordaintes for interpolation
 * @param wi w coordaintes for interpolation
 * @param method cubic or linear
 * @return interp the interpolated input matrix
 */
std::vector<arma::cube> Interpolate4D(const arma::vec& x, const arma::vec& y, const arma::vec& z, const arma::vec& w, const std::vector<arma::cube>& in, const arma::vec& xi, const arma::vec& yi, const arma::vec& zi, const arma::vec& wi, const std::string& method)
{
    //set interpolated cube
    //temporary cube for in plane interpolation
    std::vector<arma::cube> interpTemp1(w.size());
    std::vector<arma::cube> interpTemp2(w.size());
    std::vector<arma::cube> interpTemp3(w.size());
    std::vector<arma::cube> interp(wi.size());

    for(unsigned int i = 0; i < w.size(); i++)
    {
        interpTemp1[i].set_size(xi.size(), y.size(), z.size());
        interpTemp2[i].set_size(xi.size(), yi.size(), z.size());
        interpTemp3[i].set_size(xi.size(), yi.size(), zi.size());
        interp[i].set_size(xi.size(), yi.size(), zi.size());
    }

    for(unsigned int i = 0; i < wi.size(); i++)
    {
        interp[i].set_size(xi.size(), yi.size(), zi.size());
    }

    //only interpolate if necessary
    if(accu(abs(x - xi)) == 0 && accu(abs(y - yi)) == 0 && accu(abs(z - zi)) == 0)
    {
        interpTemp3 = in;
    }
    else
    {
        //interpolate in the x-y plane
        for(unsigned int i = 0; i < w.size(); i++)
        {
            for(unsigned int j = 0; j < z.size(); j++)
            {
                for(unsigned int k = 0; k < y.size(); k++)
                {
                    arma::vec initData = in[i].slice(j).col(k);

                    //holder for interpolated data
                    arma::vec interpData;

                    //interpolate data
                    //arma::interp1(z, initData, zi, interpData);
                    Interpolate1D(initData, x, interpData, xi, method);

                    //linear extrapolation if the data is outside region
                    if(interpData.has_nan() == true)
                    {
                        arma::uvec nonFiniteCoords = arma::find_nonfinite(interpData);

                        if(nonFiniteCoords[0] == 0 && nonFiniteCoords[nonFiniteCoords.size() - 1] == interpData.size() - 1)
                        {
                            interpData[0] = 2.0*interpData[1] - interpData[2];
                            interpData[interpData.size() - 1] = 2.0*interpData[interpData.size() - 2] - interpData[interpData.size() - 3];
                        }
                    }
                    else if(interpData[0] == 0 && interpData[interpData.size() - 1] == 0)
                    {
                        interpData[0] = 2.0*interpData[1] - interpData[2];
                        interpData[interpData.size() - 1] = 2.0*interpData[interpData.size() - 2] - interpData[interpData.size() - 3];
                    }

                    interpTemp1[i].slice(j).col(k) = interpData;
                }

                for(unsigned int k = 0; k < xi.size(); k++)
                {
                    arma::vec initData = interpTemp1[i].slice(j).row(k);

                    //holder for interpolated data
                    arma::vec interpData;

                    //interpolate data
                    //arma::interp1(z, initData, zi, interpData);
                    Interpolate1D(initData, y, interpData, yi, method);

                    //linear extrapolation if the data is outside region
                    if(interpData.has_nan() == true)
                    {
                        arma::uvec nonFiniteCoords = arma::find_nonfinite(interpData);

                        if(nonFiniteCoords[0] == 0 && nonFiniteCoords[nonFiniteCoords.size() - 1] == interpData.size() - 1)
                        {
                            interpData[0] = 2.0*interpData[1] - interpData[2];
                            interpData[interpData.size() - 1] = 2.0*interpData[interpData.size() - 2] - interpData[interpData.size() - 3];
                        }
                    }
                    else if(interpData[0] == 0 && interpData[interpData.size() - 1] == 0)
                    {
                        interpData[0] = 2.0*interpData[1] - interpData[2];
                        interpData[interpData.size() - 1] = 2.0*interpData[interpData.size() - 2] - interpData[interpData.size() - 3];
                    }

                    interpTemp2[i].slice(j).row(k) = interpData;
                }
            }      

            //interpolate across z axis
            for(unsigned int j = 0; j < yi.size(); j++)
            {
                for(unsigned int k = 0; k < xi.size(); k++)
                {

                    arma::vec initData = interpTemp2[i].tube(k,j);

                    //holder for interpolated data
                    arma::vec interpData;

                    //interpolate data
                    //arma::interp1(z, initData, zi, interpData);
                    Interpolate1D(initData, z, interpData, zi, method);

                    //linear extrapolation if the data is outside region
                    if(interpData.has_nan() == true)
                    {
                        arma::uvec nonFiniteCoords = arma::find_nonfinite(interpData);

                        if(nonFiniteCoords[0] == 0 && nonFiniteCoords[nonFiniteCoords.size() - 1] == interpData.size() - 1)
                        {
                            interpData[0] = 2.0*interpData[1] - interpData[2];
                            interpData[interpData.size() - 1] = 2.0*interpData[interpData.size() - 2] - interpData[interpData.size() - 3];
                        }
                    }
                    else if(interpData[0] == 0 && interpData[interpData.size() - 1] == 0)
                    {
                        interpData[0] = 2.0*interpData[1] - interpData[2];
                        interpData[interpData.size() - 1] = 2.0*interpData[interpData.size() - 2] - interpData[interpData.size() - 3];
                    }

                    interpTemp3[i].tube(k,j) = interpData;
                }
            }
        }

    }


    //interpolate across frequency 
    for(unsigned int i = 0; i < xi.size(); i++)
    {
        for(unsigned int j = 0; j < yi.size(); j++)
        {
            for(unsigned int k = 0; k < zi.size(); k++)
            {
                arma::vec initData(w.size());

                for(unsigned int l = 0; l < w.size(); l++)
                {                    
                    initData[l] = interpTemp3[l](i,j,k);
                }

                //holder for interpolated data
                arma::vec interpData;

                //interpolate data
                //arma::interp1(z, initData, zi, interpData);
                Interpolate1D(initData, w, interpData, wi, method);

                //linear extrapolation if the data is outside region
                if(interpData.has_nan() == true)
                {
                    arma::uvec nonFiniteCoords = arma::find_nonfinite(interpData);

                    if(nonFiniteCoords[0] == 0 && nonFiniteCoords[nonFiniteCoords.size() - 1] == interpData.size() - 1)
                    {
                        interpData[0] = 2.0*interpData[1] - interpData[2];
                        interpData[interpData.size() - 1] = 2.0*interpData[interpData.size() - 2] - interpData[interpData.size() - 3];
                    }
                }
                else if(interpData[0] == 0 && interpData[interpData.size() - 1] == 0)
                {
                    interpData[0] = 2.0*interpData[1] - interpData[2];
                    interpData[interpData.size() - 1] = 2.0*interpData[interpData.size() - 2] - interpData[interpData.size() - 3];
                }

                for(unsigned int l = 0; l < wi.size(); l++)
                {  
                    interp[l](i,j,k) = interpData[l];
                }
            }
        }
    }

    return interp;

}

/**
 * @brief Function to flatten a vector of cubes into one vecotr
 * 
 * @param in a vector of cubes
 * @return out: a column vector
 */
arma::cx_vec Flatten(const std::vector<arma::cx_cube>& in)
{
    //get the number of elemetns in the input
    int nElemSlice = in[0].n_elem;
    int nElem = in.size()*nElemSlice;

    //output vec
    arma::cx_vec out(nElem);

    for(unsigned int i = 0; i < in.size(); i++)
    {
        arma::cx_vec temp = arma::vectorise(in[i]);
        out.subvec(i*nElemSlice, size(temp)) = temp;
    }

    return out;
}

/**
 * @brief Function to make a column vector 4D, matching the size of another 4D object
 * 
 * @param in a column vector
 * @param sizeMatch system to match the size d 
 * @return out: 4d matrix
 */
std::vector<arma::cx_cube> Make4D(const arma::cx_vec& in, const std::vector<arma::cx_cube>& sizeMatch)
{
    //set size of output
    int length = sizeMatch.size();
    std::vector<arma::cx_cube> out(length);
    int x = sizeMatch[0].n_cols;
    int y = sizeMatch[0].n_rows;
    int z = sizeMatch[0].n_slices;

    //allocate values
    for(int i = 0; i < length; i++)
    {
        out[i].set_size(x, y, z);
        for(int j = 0; j < z; j++)
        {
            for(int k = 0; k < y; k++)
            {
                for(int l = 0; l < x; l++)
                {
                    out[i](l, k, j) = in[l + x*k + x*y*j + x*y*x*i];
                }
            }
        }
    }

    return out;
}

/**
 * @brief General permutation function for a cube
 * 
 * @param in cube to be permuted
 * @param order order of permutation
 * @return out: permuted cube
 */
arma::cx_cube Permute(const arma::cx_cube& in, const std::string& order)
{
    arma::cx_cube out;
    arma::cx_cube tmp;

    if(order == "132")
    {
        out = permute23Comp(in);
    }
    else if(order == "213")
    {
        //regular tranpose
        //tranpose each slice within the cube
        out = in;
        out.each_slice([](arma::cx_mat& tempA){tempA = tempA.st();});
    }
    else if(order == "231")
    {
        //regular tranpose + permute 23
        //tranpose each slice within the cube
        tmp = in;
        tmp.each_slice([](arma::cx_mat& tempA){tempA = tempA.st();});
        out = permute23Comp(tmp);
    }
    else if(order == "321")
    {
        //regular tranpose + permute 23
        //tranpose each slice within the cube
        tmp = in;
        tmp.each_slice([](arma::cx_mat& tempA){tempA = tempA.st();});
        out = permute23Comp(tmp);
        out.each_slice([](arma::cx_mat& tempA){tempA = tempA.st();});
    }
    else if(order == "312")
    {
        //regular tranpose + permute 23
        //tranpose each slice within the cube
        tmp = in;
        out = permute23Comp(tmp);
        out.each_slice([](arma::cx_mat& tempA){tempA = tempA.st();});
    }
    else
    {
        std::cout << "Undefined permutation order" << std::endl;
        exit(1);
    }

    return out;
}


/**
 * @brief Template to determine the sign of a number
 * 
 * @tparam T 
 * @param x
 * @return int 
 */
template <typename T> int sgn(T x) 
{
    return (T(0) < x) - (x < T(0));
}


/**
 * @brief Permutation template for cubes
 * 
 * @tparam T 
 * @param cube 
 * @return Cube<T> 
 */
template <typename T>
static arma::Cube<T> permute23(const arma::Cube<T> &cube)
{
    const arma::uword d1 = cube.n_rows;
    const arma::uword d2 = cube.n_cols;
    const arma::uword d3 = cube.n_slices;
    const arma::uword d1TimesD3 = d1 * d3;
    const arma::Cube<T> output(d1, d3, d2);

    const T *from = cube.memptr();
    T *to = (double*)(output.memptr());

    for (arma::uword s = 0; s < d3; ++s){
        T *tmp = to + d1 * s;
        for (arma::uword c = 0; c < d2; ++c){
            memcpy(tmp, from, d1 * sizeof(*from));
            from += d1;
            tmp += d1TimesD3;
        }
    }

    return output;
}


/**
 * @brief Permutation template for cubes
 * 
 * @tparam T 
 * @param cube 
 * @return Cube<T> 
 */
template <typename T>
static arma::Cube<T> permute23Comp(const arma::Cube<T> &cube)
{
    const arma::uword d1 = cube.n_rows;
    const arma::uword d2 = cube.n_cols;
    const arma::uword d3 = cube.n_slices;
    const arma::uword d1TimesD3 = d1 * d3;
    const arma::Cube<T> output(d1, d3, d2);

    const T *from = cube.memptr();
    T *to = (std::complex<double>*)(output.memptr());

    for (arma::uword s = 0; s < d3; ++s){
        T *tmp = to + d1 * s;
        for (arma::uword c = 0; c < d2; ++c){
            memcpy(tmp, from, d1 * sizeof(*from));
            from += d1;
            tmp += d1TimesD3;
        }
    }

    return output;
}


/**
 * @brief Extern capabilities to allow interfacing with python code
 * 
 */
extern "C" {
    Eliashberg* Eliashberg_New(const char* config)
    {
        //convert char* to string
        std::string strConfig(config);   
        return new Eliashberg(strConfig); 

    }; 

    void Eliashberg_SolveEliashberg(Eliashberg* eliashberg){ eliashberg->SolveEliashberg(); }
}





