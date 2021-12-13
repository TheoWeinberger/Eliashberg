#--------------------------------------------------------------------------
# Spectrum.py - Code to fit a functional form of the spin fluctuation
# susceptibility uses least-squares regression and return the fits to 
# of the parameters, this data can then be used by the Eliashberg
# code to run an Eliashberg calculation
#
# Ideally should do a global search in the whole of q,w space to get the 
# parameters that fit the whole system
# However, if that global fit fails, a fit of the static susceptibility
# is performed on the central w data and then the rest of the parameters
# fitten in omega space
#
# For the input, the number of peaks and their location should be 
# specified for the moment, in the future this should be automised
# consider using moving average/wavelets for this 
# https://mathematica.stackexchange.com/questions/23828/finding-local-minima-maxima-in-noisy-data
#
# version: 1
# date: 02/12/2021
# author: Theo Weinberger
#--------------------------------------------------------------------------

#import scipy.optimize to use curve_fit for optimsation
from os import error
import scipy.optimize as optimize
#lmfit might be better
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#see
# https://stackoverflow.com/questions/12208634/fitting-only-one-parameter-of-a-function-with-many-parameters-in-python
# https://stackoverflow.com/questions/34136737/using-scipy-curve-fit-for-a-variable-number-of-parameters

def omega_func(omega, A, chi_q_inv, gamma_q, omega_0):
    """
    This function is used to find the variable parameters
    chi_q_inv and gamma_q. It returns chi_imag which is the imaginary part
    of the spin response function. It is used in the case that suscpetibility 
    measurements are taken only as a function of energy/frequency i.e. at just
    one q vector. Therefore this will give chi_q_inv gamma_q which can subsequently
    be evaluated using the relationship chi_q_inv = chi_Q_inv + c|q-Q|^2 where 
    if measurements are taken at q multiple points at a given frequency.

    For AFM gamma_q = gamma while for FM gamma_q = q.gamma

    Args:

        omega: the frequency of incident neutrons

        A: normalisation

        chi_q_inv: the inverse of the static spin susceptibility

        gamma_q: the q-dependent damping

        omega_0: location of frequency peak

    Returns:

        chi_int: the neutron intensity
    """

    chi_int = -A*(omega*gamma_q)/((gamma_q*chi_q_inv)**2 + (omega-omega_0)**2)

    return chi_int


def chi_imag_AFM_func(coords, chi_Q_inv, gamma, c, Q, freq_0):
    """
    This function is used to calculate the antiferromagnetic contribution 
    to the imaginary part of the response function

    Args:

        coords: an 2d array containing the q space coords and the frequency 
        coords order is [q,w]

        chi_Q_inv: the maximal value of the static susceptibility peak

        gamma: the damping constant

        c: the decay factor in the static susceptibility

        Q: the Q coordinate of the susceptibility peak 

        freq_0: location of frequency peak

    Returns:
    
        chi_imag: the imaginary part of the response function
    """

    #extract frequency and momentum components
    freq = coords[:,1]
    q = coords[:,0]

    chi_imag = (freq*gamma)/((gamma*(chi_Q_inv + c*(q-Q)**2))**2 + (freq - freq_0)**2)

    return chi_imag


def chi_imag_FM_func(coords, chi_Q_inv, gamma, c, Q, freq_0):
    """
    This function is used to calculate the ferromagnetic contribution 
    to the imaginary part of the response function

    Args:

        coords: an 2d array containing the q space coords and the frequency 
        coords order is [q,w]

        chi_Q_inv: the maximal value of the static susceptibility peak

        gamma: the damping constant

        c: the decay factor in the static susceptibility

        Q: the Q coordinate of the susceptibility peak 

        freq_0: location of frequency peak

    Returns:
    
        chi_imag: the imaginary part of the response function
    """

    #extract frequency and momentum components
    freq = coords[:,1]
    q = coords[:,0]

    temp1 = (freq*q*gamma)
    temp2 = ((gamma*q*(chi_Q_inv + c*(q-Q)**2))**2 + (freq-freq_0)**2)

    chi_imag = np.divide(temp1, temp2, out=np.zeros_like(temp1), where=temp2!=0.0)

    return chi_imag


def intensity_func(coords, A, chi_Q_inv, gamma, c, Q, freq_0, mag_model, num_peaks):
    """
    This functions sums over the responses for the individual peaks
    to calculate the overall neutron intensity spectrum as a function 
    of frequency and q

    Args:

        coords: an 2d array containing the q space coords and the frequency 
        coords order is [q,w]

        A: a scaling parameter

        chi_Q_inv: the maximal value of the static susceptibility peak

        gamma: the damping constant

        c: the decay factor in the static susceptibility

        Q: the Q coordinate of the susceptibility peak 

        freq_0: location of frequency peak

        mag_model: the type of magnetism the peak corresponds to

        num_peaks: the number of peaks in the system

    Return:

        intesity: the neutron intensity

    """
    #enforce positivity
    chi_Q_inv = [abs(val) for val in chi_Q_inv]
    gamma = [abs(val) for val in gamma]
    c = [abs(val) for val in c]

    intensity = 0

    for i in range(num_peaks):

        if(mag_model[i] == 'FM'):

            intensity += -A*chi_imag_FM_func(coords, chi_Q_inv[i], gamma[i], c[i], Q[i], freq_0[i])

        elif(mag_model[i] == 'AFM'):

            intensity += -A*chi_imag_AFM_func(coords, chi_Q_inv[i], gamma[i], c[i], Q[i], freq_0[i])

        else:

            print("Wrong inputs in the magnetic model list, please select either FM or AFM")
            exit()

    return intensity


def wrapper_intensity_func(coords, mag_model, num_peaks, Q, *args):
    """
    A wrapper function to convert the input for the intensity function 
    data from a list of all params in args to a lists for each parameter

    Args:

        coords: an 2d array containing the q space coords and the frequency 
        coords order is [q,w]

        mag_model: the type of magnetism the peak corresponds to

        num_peaks: the number of peaks in the system

        Q: the peak positions in Q-space

        args: system parameters in single list format

    Return:

        intesity: the neutron intensity

    """
    #extract variable parameters
    A = args[0][0]
    chi_Q_inv, gamma, c, freq_0 = list(args[0][1:num_peaks + 1]), list(args[0][num_peaks + 1: 2*num_peaks + 1]), list(args[0][2*num_peaks + 1: 3*num_peaks + 1]), list(args[0][3*num_peaks + 1: 4*num_peaks + 1])

    #calculate intensity
    intensity = intensity_func(coords, A, chi_Q_inv, gamma, c, Q, freq_0, mag_model, num_peaks)

    return intensity


#main script
if __name__ == "__main__":

    #read in data
    #read in raw data
    df = pd.read_csv("YFe2Ge2.dat", sep = ",", engine="python")

    #extract input data
    data = df["Intensity"].to_numpy()
    errors = df["Error"].to_numpy()
    #protect against divide by zero
    errors = np.divide(errors, abs(data), out=np.zeros_like(errors), where=data>=0.1)
    errors[errors == 0] = np.partition(np.unique(errors),1)[1]
    
    q = df["Q"].to_numpy()
    e = df["Energy"].to_numpy()

    #Data params
    num_peaks = 2
    #peak types
    mag_model = ['AFM', 'AFM']
    #peak positions
    Q = [0.0, 0.5]

    #initial params
    A = [-2.0]
    chi_Q_inv = [1.0, 1.0]
    gamma = [6.0, 1.0]
    c = [100.0, 100.0]
    freq_0 = [1.0, 1.0]

    #create parameters list
    params = A + chi_Q_inv + gamma + c + freq_0

    #reformat coords for fitting
    coords = np.array([q,e], dtype=object)
    coords = coords.T

    #fit curve
    popt, pcov = optimize.curve_fit(lambda coords, *params_init: wrapper_intensity_func(coords, mag_model, num_peaks, Q, params_init), coords, data, p0=params)

    #reformat fitted parameters
    popt = abs(popt)
    popt[0] = -popt[0]

    #output parameters
    print(popt)

    #plot fitted results at on energy value
    e = np.full((1000),6)
    q = np.linspace(-0.25,0.75,1000)

    coords = np.array([q,e], dtype=object)
    coords = coords.T
    intensity = wrapper_intensity_func(coords, mag_model, num_peaks, Q, popt)

    plt.plot(q, intensity)
    plt.errorbar(df["Q"], df["Intensity"], df["Error"])
    plt.show()

    """
    #fit in just energy space and compare the fits
    #read in data
    #read in raw data
    df = pd.read_csv("YFe2Ge2EQ0_Clean.dat", sep = ",", engine="python")

    #extract input data
    data = df["Intensity"].to_numpy()
    errors = df["Error"].to_numpy()
    #protect against divide by zero
    errors = np.divide(errors, abs(data), out=np.zeros_like(errors), where=data>=0.1)
    errors[errors == 0] = np.partition(np.unique(errors),1)[1]
    
    q = df["Q"].to_numpy()
    e = df["Energy"].to_numpy()

    #Data params
    num_peaks = 1
    #peak types
    mag_model = ['AFM']
    #peak positions
    Q = [0.0]

    #initial params
    A = [-2.0]
    chi_q_inv = [1.0]
    gamma_q = [6.0]

    #create parameters list
    params = A + chi_q_inv + gamma_q + gamma_q

    #reformat coords for fitting
    coords = np.array([q,e], dtype=object)
    coords = coords.T

    #fit curve
    popt, pcov = optimize.curve_fit(omega_func, coords[:,1], data, p0=params, maxfev = 1000000)

    #reformat fitted parameters
    popt = abs(popt)
    popt[0] = -popt[0]

    #output parameters
    print(popt)

    #plot fitted results at on energy value
    q = np.full((1000),0.0)
    e = np.linspace(0,10,1000)

    coords = np.array([q,e], dtype=object)
    coords = coords.T
    intensity = omega_func(coords[:,1], *popt)

    plt.plot(e, intensity)
    plt.errorbar(df["Energy"], df["Intensity"], df["Error"])
    plt.show()
    """














    








