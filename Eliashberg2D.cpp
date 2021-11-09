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


#include "Eliashberg2D.hpp"
#include "Eliashberg2DSettings.hpp"


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

    //set default plotting value
    _plot = "g";


    std::cout << " " << std::endl;

    ReadFile(fileName, _t, _ratioTight, _a, _doping, _magModel, _phiModel, _chi0,
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

    int lenKSample = _kSquaredSample.size();
    int lenGSample = _gSquaredChi0tSample.size();

    //check sampling is properly set up
    if(lenKSample > 1 && _plot != "k")
    {
        std::cout << "Data structure is not formatted for plotting in k space, please check the k sampling matrix and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of k, _plot should be set to k and the k sampling array should have multiple values" << std::endl;
        exit(1);
    }

    if(lenKSample > 1 && _plot != "k")
    {
        std::cout << "Data structure is not formatted for plotting in g space, please check the g sampling matrix and the plotting settings" << std::endl;
        std::cout << "To plot Tc as a function of g, _plot should be set to g and the g sampling array should have multiple values" << std::endl;
        exit(1);
    }

    //temperature arrays
    arma::mat tC(lenKSample, lenGSample, arma::fill::zeros);
    /*********************************/
    arma::mat tInit(lenKSample, lenGSample, arma::fill::zeros); //does this get used
    /**********************************/

    //calculate the array containing the chemical potential
    arma::vec muArray = _calcMu(_t0); //note 0 is picked as an initial trial value out of convenience

    //direct test against Ran code
    //muArray = {1.351850567693204e+03, 6.765567525297223e+02, 3.405198980389718e+02, 1.765050769606098e+02, 1.033043515584506e+02, 86.539442523009030, 1.030800627209959e+02, 1.274838994943234e+02, 1.452422053806421e+02, 1.541559644533966e+02, 1.582541906946000e+02, 1.596428565534779e+02, 1.604206942145785e+02};

    for(int a = 0; a < lenKSample; a++)
    {
        for(int b = 0; b < lenGSample; b++)
        {
            
            //user information
            std::cout << "Initialising run with parameters" << std::endl;
            std::cout << "g: " << _gSquaredChi0tSample[b] << std::endl;
            std::cout << "k: " << _kSquaredSample[a] << std::endl;

            //initialise energy
            _energy = _Dispersion(qX, qY);


            //initialise parameters for search
            _kSquared = _kSquaredSample[a];
            double gSquaredChi0t = _gSquaredChi0tSample[b];
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
            else
            {
                std::cout << "Invalid value of the magnetic model" << std::endl;
                std::cout << "Please enter either FM or AFM" << std::endl;
                exit(1);
            }

            /*initialise the relaxation weight index
            The relaxation weight smooths the transition between 
            steps to allow for more stable convergence to self-consistency*/
            unsigned int relaxIndex = 0;

            /*Vectors containing the lambda values and temperature steps
            these start off as being unit length and are appended to each 
            iterations. Note the appending operation is highly inefficient 
            and so should be altered in later code*/
            arma::vec tStep(1);
            arma::vec lambdaVec(1);

            //doubles for temperature and sampling number
            double t = _t0;
            double n = _n0;

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
                t = _t0;
                n = _n0;

                unsigned int muIndex = 0;

                //allocate memory for sigma and phi
                arma::cx_cube sigma(_nK, _nK, 2*_n0, arma::fill::zeros);
                arma::cube phi(_nK, _nK, 2*_n0);
                arma::cube phiFilter(_nK, _nK, 2*_n0);

                //initialise the RG corrections 
                arma::cx_cube dSigma(_nK, _nK, 2*_n0, arma::fill::zeros); 
                arma::cube dPhi(_nK, _nK, 2*_n0, arma::fill::zeros); 

                //counter for iteration number
                int counter = 0;

                //boolean to see if it is the first run
                bool firstRun;

                //boolean to state whether regular matsubara convolution is occuring or the RG corrections
                bool rg = false;

                /*Search for self consistency  of the Green's function
                is for when the value of lambda tends to 1. Therefore
                the search ends when lambda = 1*/
                while(abs(lambdaVec(lambdaVec.size() - 1)) < 1.0)
                {

                    //reset all cubes to be filled with 0
                    sigma.zeros();
                    phi.zeros();
                    phiFilter.zeros();

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

                    for(unsigned int i = 0; i < qX.n_elem; i++)
                    {
                        for(unsigned int j = 0; j < qY.n_elem; j++)
                        {
                            for(unsigned int k = 0; k < wMatsu.n_elem; k++)
                            {
                                gInv(i,j,k) = (I*wMatsu(k) - (_energy(i,j) - _mu));
                            }
                        }
                    }

                    g = 1/gInv;

                    //Calculate chinu data, this corresponds to the analystical intergral of the dynamical suscpetibility
                    arma::cube chiQNu =  _ChiQNu(vMatsu, qX, qY);
                    arma::cx_cube chiQNuComplex = RealToComplex(chiQNu);

                    //double to store the relative error in the convergence in sigma
                    double relErrS = INFINITY;

                    //solve for sigma
                    for(int i = 0; i < _maxIter; i++)
                    {
                        //evaluate new Greens function
                        g = 1/(gInv - sigma);

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
                        arma::cx_cube sigmaMatsuConv = gSquared*t/(pow(_nK, 2.0))*_MatsuConv(chiQNuComplex, g, 2*n, 4*n - 1, firstRun, rg) + dSigma;

                        //set firstRun to false for remaining iterations
                        firstRun = false;

                        //work out relative erorr in convergence
                        double deltaS = arma::accu(abs(sigmaMatsuConv-sigma));
                        double totalS = arma::accu(abs(sigma));
                        relErrS = deltaS/totalS;

                        //iterate the sigma matrix
                        sigma = _relaxation[relaxIndex]*sigmaMatsuConv + (1 - _relaxation[relaxIndex])*sigma;                     
                    
                        /*sigma should be symmetric in rows and columns (i.e. within the slices)
                        make sure that this is the case as it can vary a bit over time due to 
                        floating point/rounding errors*/
                        if(_symm == true)
                        {
                            sigma = Symmetrise(sigma);
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
                    double scale = abs(dPhi).max();

                    //scale phi if needed
                    if(scale > 0)
                    {
                        phi *= scale*_phiRatio;
                    }

                    //calcualte <Phi|Phi> , shoudln;t this be a proper inner product???
                    /*************************************************************/
                    double phiInner = arma::accu(phi%phi);

                    //get the absolute value of g
                    arma::cx_cube gAbsSquaredComp = RealToComplex(pow(abs(g), 2.0));

                    //convolve chi wth |G|^2%phi
                    arma::cx_cube gAbsPhi = gAbsSquaredComp%phi;

                    /********************
                     * 
                     * in general phimatsu conv can be imaginary, change this
                     **********************************/
                    arma::cube phiMatsuConv = arma::real(coupling*t/pow(_nK, 2.0)*_MatsuConv(chiQNuComplex, gAbsPhi, 2*n, 4*n - 1, firstRun, rg)) + dPhi;

                    //calcilate lambda as <Phi|Phi> = <Phi|A|Phi> = <Phi|Phi1>
                    double lambda = arma::accu(phi%phiMatsuConv)/phiInner;

                    //apply symmetry transforms
                    phi = _SymmByFiltLabel(phiMatsuConv, phiFilter);

                    //initialise relative error and normalisation for lambda search
                    double relErrL = INFINITY;
                    bool normalise = false;

                    for(int i = 0; i < _maxIter; i++)
                    {

                        //calcualte <Phi|Phi>
                        phiInner = arma::accu(phi%phi);


                        //convolve chi wth |G|^2%phi
                        gAbsPhi = gAbsSquaredComp%phi;
                        phiMatsuConv = arma::real(coupling*t/pow(_nK, 2.0)*_MatsuConv(chiQNuComplex, gAbsPhi, 2*n, 4*n - 1, firstRun, rg)) + dPhi;

                        //calcilate lambda as <Phi|Phi> = <Phi|A|Phi> = <Phi|Phi1>
                        double lambda1 = arma::accu(abs(phi%phiMatsuConv))/phiInner;

                        //calculate the relative change in lambda
                        double relErrL1 = abs(lambda1 - lambda)/abs(lambda);
                        
                        //check normalisation cases to update lambda
                        if(normalise == true)
                        {
                            relErrL = relErrL1;

                            //iterate to the next step
                            phi = _SymmByFiltLabel(phiMatsuConv, phiFilter)/lambda1;
                            lambda = lambda1;
                        }
                        else if(relErrL1 > relErrL)
                        {
                            //if relative error increases, set to normalise
                            normalise = true;
                            phi = phi/lambda;
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
                        phi /= lambda;
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
                    arma::cube chiCut =  _ChiQNu(vCut, qX, qY);
                    arma::cx_cube chiCutComplex = RealToComplex(chiCut);

                    arma::cx_cube gL = g.slices((n - m),(n + m -1));

                    rg = true;

                    //calculate the contribution to this iteration from L to L domains
                    arma::cx_cube sigmaLL = gSquared*t/(pow(_nK, 2.0))*_MatsuConv(chiCutComplex, gL, 2*m, 4*m - 1, firstRun, rg);

                    //contribution to sigma from L to H domains
                    arma::cx_cube sigmaL = sigma.slices((n - m),(n + m -1)); 

                    //holder for dSigma
                    arma::cube dSigmaReal, dSigmaImag;

                    //set to zero
                    dPhi.zeros();
                    dSigma.zeros();
                    dSigmaImag.zeros();
                    dSigmaReal.zeros();


                    //note I think 3D interpolation is unecessary as it is just along the frequency domain that interpolation occurs
                    //interpolate real and imaginary parts
                    dSigmaReal = Interpolate3D(qX, qY, wL, arma::real(sigmaL - sigmaLL), qX, qY, wMatsu/2.0);
                    dSigmaImag = Interpolate3D(qX, qY, wL, arma::imag(sigmaL - sigmaLL), qX, qY, wMatsu/2.0);

                    //set whole dSigma matrix
                    dSigma.copy_size(dSigmaReal);
                    dSigma.set_real(dSigmaReal);
                    dSigma.set_imag(dSigmaImag);

                    //perform same cutting for phi

                    arma::cube phiL = phi.slices((n - m),(n + m - 1)); 

                    arma::cx_cube gAbsPhiL = RealToComplex(pow(abs(gL), 2.0))%phiL;

                    arma::cube phiLL = arma::real(coupling*t/pow(_nK, 2.0)*_MatsuConv(chiCutComplex, gAbsPhiL, 2*m, 4*m - 1, firstRun, rg));

                    dPhi = Interpolate3D(qX, qY, wL, lambda*phiL - phiLL, qX, qY, wMatsu/2.0);

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
            tC(a, b) = tInterp(indexTC);
            //output the critical temperature
            std::cout << "    " << std::endl;
            std::cout << "Tc: "<< tC(a,b) << std::endl;
            std::cout << "Tc/Tsf: " << tC(a,b)/_tSF << std::endl;
            std::cout << "    " << std::endl;

        }      
    } 

    //output k/g data if relevant
    if(lenKSample > 1)
    {
        std::cout << "k: " << std::endl;
        _kSquaredSample.t().print();

        //output all the critical temperature
        std::cout << "Tc: " << std::endl;
        tC.t().print();
        std::cout << "Tc/Tsf: " << std::endl;
        (tC/_tSF).t().print();
    }

    if(lenGSample > 1)
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

        outputData.save("gData", arma::csv_ascii);
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

        outputData.save("kData", arma::csv_ascii);
    }
    else
    {
        std::cout << "Incorrect plotting values, please set plot to either g or k" << std::endl;
        exit(1);
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

    int iMax = qX.n_elem;
    int jMax = qY.n_elem;
    int kMax = vMatsu.n_elem;

    arma::cube chiQNu(iMax, jMax, kMax);

    double omega0, eta, qHat, a;

    //Evaluate the integral of chiQNu at each point
    for(int i = 0; i < iMax; i++)
    {
        for(int j = 0; j < jMax; j++)
        {
            for(int k = 0; k < kMax; k++)
            {

                omega0 = _Omega0(qX[i], qY[j]);
                eta = _Eta(qX[i], qY[j]);

                //qhat value depends on magnetic model
                if(_magModel == "FM")
                {
                    qHat = _QMinus(qX[i], qY[j]);
                }
                else if(_magModel == "AFM")
                {
                    qHat = _QPlus(qX[i], qY[j]);
                }
                else
                {
                    std::cout << "Invalid value of the magnetic model" << std::endl;
                    std::cout << "Please enter either FM or AFM" << std::endl;
                    exit(1);
                }

                //last term in the integral of the dynamics susceptibility
                a = eta*(_kSquared + pow(qHat, 2.0));

                chiQNu(i, j, k) = (2*_chi0/M_PI)*_Omega0(qX[i], qY[j])*_ChiInt(omega0, vMatsu[k], a);

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
 * @param firstRun boolean specifying if it is the first
 * @param boolean specifying if rg corrections are occuring 
 * @return matrixConv The matrix corresponding to a linear convolution in frequency of A and B
 */
arma::cx_cube Eliashberg::_MatsuConv(const arma::cx_cube& matrixA, const arma::cx_cube& matrixB, const int& lowIndex, const int& highIndex, const bool& firstRun, const bool& rg)
{

    //These tempoarary cubes will store the padded data 
    arma::cx_cube matrixAPadded = matrixA;
    arma::cx_cube matrixBPadded = matrixB;

    //Apply circular shifts to get frequency centred output
    matrixAPadded = IfftShift(matrixAPadded, 1);
    matrixAPadded = IfftShift(matrixAPadded, 2);

    matrixBPadded = IfftShift(matrixBPadded, 1);
    matrixBPadded = IfftShift(matrixBPadded, 2); 

    //pad matrices
    matrixAPadded = PadCube(matrixAPadded, matrixB.n_slices - 1, 0.0);
    matrixBPadded = PadCube(matrixBPadded, matrixA.n_slices - 1, 0.0);

    arma::cx_cube convolution;

    if(rg == true)
    {


        //cube for the planner to use and overwrite
        arma::cx_cube plannerCube;

        plannerCube.copy_size(matrixAPadded);

        _SetDFTPlansRG(plannerCube, plannerCube);

        //apply ffts, multiple and then apply ifft to achieve convolution
        fftw_execute_dft(_forwardPlanRG, (double(*)[2])&matrixAPadded(0,0,0), (double(*)[2])&matrixAPadded(0,0,0));
        fftw_execute_dft(_forwardPlanRG, (double(*)[2])&matrixBPadded(0,0,0), (double(*)[2])&matrixBPadded(0,0,0));
        

        convolution = matrixAPadded%matrixBPadded;

        fftw_execute_dft(_inversePlanRG, (double(*)[2])&convolution(0,0,0), (double(*)[2])&convolution(0,0,0));
    }
    else
    {

        //set DFTs if it is the first run
        if(firstRun == true)
        {
            //cube for the planner to use and overwrite
            arma::cx_cube plannerCube;

            plannerCube.copy_size(matrixAPadded);

            _SetDFTPlans(plannerCube, plannerCube);
        }

        //apply ffts, multiple and then apply ifft to achieve convolution
        fftw_execute_dft(_forwardPlan, (double(*)[2])&matrixAPadded(0,0,0), (double(*)[2])&matrixAPadded(0,0,0));
        fftw_execute_dft(_forwardPlan, (double(*)[2])&matrixBPadded(0,0,0), (double(*)[2])&matrixBPadded(0,0,0));

        convolution = matrixAPadded%matrixBPadded;

        fftw_execute_dft(_inversePlan, (double(*)[2])&convolution(0,0,0), (double(*)[2])&convolution(0,0,0));

    }

    convolution /= convolution.n_elem;

    //truncate the frequency domain to remain within the cutoff range
    arma::cx_cube convolutionCrop = convolution.slices(lowIndex - 1, highIndex - 1);

    //reshift matrix
    convolutionCrop = IfftShift(convolutionCrop, 1);
    convolutionCrop = IfftShift(convolutionCrop, 2);


    return convolutionCrop; 

}


/**
 * @brief Function to set DFT plans for the matsubara frequency convolution
 * 
 * @param in The matrix being transformed
 * @param out The output matrix of the DFT
 */
void Eliashberg::_SetDFTPlans(const arma::cx_cube& in, const arma::cx_cube& out)
{

    _forwardPlan = fftw_plan_dft_3d(in.n_slices, in.n_cols, in.n_rows, (double(*)[2])&in(0,0,0), (double(*)[2])&out(0,0,0), FFTW_FORWARD, FFTW_MEASURE);

    _inversePlan = fftw_plan_dft_3d(in.n_slices, in.n_cols, in.n_rows, (double(*)[2])&in(0,0,0), (double(*)[2])&out(0,0,0), FFTW_BACKWARD, FFTW_MEASURE);

}


/**
 * @brief Function to set DFT plans for the matsubara frequency convolution
 * 
 * @param in The matrix being transformed
 * @param out The output matrix of the DFT
 */
void Eliashberg::_SetDFTPlansRG(const arma::cx_cube& in, const arma::cx_cube& out)
{

    _forwardPlanRG = fftw_plan_dft_3d(in.n_slices, in.n_cols, in.n_rows, (double(*)[2])&in(0,0,0), (double(*)[2])&out(0,0,0), FFTW_FORWARD, FFTW_ESTIMATE);

    _inversePlanRG = fftw_plan_dft_3d(in.n_slices, in.n_cols, in.n_rows, (double(*)[2])&in(0,0,0), (double(*)[2])&out(0,0,0), FFTW_BACKWARD, FFTW_ESTIMATE);

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
arma::cube Eliashberg::_PhiFun(const arma::vec& qX, const arma::vec& qY)
{

    arma::cube phi(_nK, _nK, 2*_n0);

    //case for s wave
    if(_phiModel == "s")
    {
        phi.fill(1.0);
    }
    //case for p wave
    else if(_phiModel == "p")
    {
        for(unsigned int i = 0; i < phi.n_slices; i++)
        {
            for(unsigned int j = 0; j < phi.n_cols; j++)
            {
                for(unsigned int k = 0; k < phi.n_rows; k++)
                {

                    phi(k, j, i) = sin(qX[k]*_a);
                    
                }
            }
        }        
    }
    //case for d wave
    else if(_phiModel == "d")
    {

        for(unsigned int i = 0; i < phi.n_slices; i++)
        {
            for(unsigned int j = 0; j < phi.n_cols; j++)
            {
                for(unsigned int k = 0; k < phi.n_rows; k++)
                {

                    phi(k, j, i) = cos(qX[k]*_a) -  cos(qY[j]*_a);
                    
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
arma::cube Eliashberg::_PhiSymm(const arma::vec& qX, const arma::vec& qY)
{

    arma::cube phi(_nK, _nK, 2*_n0);

    //case for s wave
    if(_phiModel == "s")
    {
        phi.fill(1.0);
    }
    //case for p wave
    else if(_phiModel == "p")
    {
        for(unsigned int i = 0; i < phi.n_slices; i++)
        {
            for(unsigned int j = 0; j < phi.n_cols; j++)
            {
                for(unsigned int k = 0; k < phi.n_rows; k++)
                {

                    phi(k, j, i) = sgn(qX[k]);
                    
                }
            }
        }        
    }
    //case for d wave
    else if(_phiModel == "d")
    {

        for(unsigned int i = 0; i < phi.n_slices; i++)
        {
            for(unsigned int j = 0; j < phi.n_cols; j++)
            {
                for(unsigned int k = 0; k < phi.n_rows; k++)
                {

                    phi(k, j, i) = sgn(qX[k] - qY[j])*sgn(qX[k] + qY[j]);
                    
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
arma::cube Eliashberg::_SymmByFiltLabel(arma::cube& matrixA, const arma::cube& filter)
{

    //initialise symmetrised matrix
    arma::cube matrixS;

    //apply symmetry transforms depending on phimodel
    if(_phiModel == "s")
    {
        matrixS = matrixA;
    }
    else if(_phiModel == "p")
    {
        arma::cube matrixA1 = matrixA%filter;
        arma::cube matrixA2 = FlipUDCube(matrixA1);
        matrixS = (matrixA1 + matrixA2)%filter/2.0;
    }
    else if(_phiModel == "d")
    {
        arma::cube matrixA1 = matrixA%filter;
        arma::cube matrixA2 = FlipLRCube(matrixA1);
        arma::cube matrixA3 = FlipUDCube(matrixA1);
        arma::cube matrixA4 = FlipLRCube(matrixA3);
        arma::cube matrixA5 = Transpose(matrixA1);
        arma::cube matrixA6 = FlipLRCube(matrixA5);
        arma::cube matrixA7 = FlipUDCube(matrixA5);
        arma::cube matrixA8 = FlipLRCube(matrixA7);
        matrixS = (matrixA1 + matrixA2 + matrixA3 + matrixA4 + matrixA5 + matrixA6 + matrixA7 + matrixA8)%filter/8.0;
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
        arma::cx_cube tempA1 = a.slices(0, shiftMax);
        arma::cx_cube tempA2 = a.slices(shiftMax + 1, length - 1);

        //create new cube with swapped slices
        arma::cx_cube shiftATemp =  arma::join_slices(tempA1, tempA2);

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
 * @return interp the interpolated input matrix
 */
arma::cube Interpolate3D(const arma::vec& x, const arma::vec& y, const arma::vec& z, const arma::cube& in, const arma::vec& xi, const arma::vec& yi, const arma::vec& zi)
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
            //arma::interp1(z, initData, zi, interpData);
            Interpolate1D(initData, z, interpData, zi, "cubic");

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




