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
import scipy.optimize as optimize
#lmfit might be better
import numpy as np
import matplotlib.pyplot as plt

#see
# https://stackoverflow.com/questions/12208634/fitting-only-one-parameter-of-a-function-with-many-parameters-in-python
# https://stackoverflow.com/questions/34136737/using-scipy-curve-fit-for-a-variable-number-of-parameters

def omega_func(omega, chi_q_inv, gamma_q):
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

        chi_q_inv: the inverse of the static spin susceptibility

        gamma_q: the q-dependent damping

    Returns:

        chi_imag: the imaginary part of the response function
    """

    chi_imag = -(omega*gamma_q)/((gamma_q*chi_q_inv)**2 + omega**2)

    return chi_imag


def chi_imag_AFM_func(coords, chi_Q_inv, gamma, c, Q):
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

    Returns:
    
        chi_imag: the imaginary part of the response function
    """

    #extract frequency and momentum components
    freq = coords[:,1]
    q = coords[:,0]

    chi_imag = (freq*gamma)/((gamma*(chi_Q_inv + c*(q-Q)**2))**2 + freq**2)

    return chi_imag


def chi_imag_FM_func(coords, chi_Q_inv, gamma, c, Q):
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

    Returns:
    
        chi_imag: the imaginary part of the response function
    """

    #extract frequency and momentum components
    freq = coords[:,1]
    q = coords[:,0]

    chi_imag = (freq*q*gamma)/((gamma*q*(chi_Q_inv + c*(q-Q)**2))**2 + freq**2)

    return chi_imag


def intensity_func(coords, A, chi_Q_inv, gamma, c, Q, mag_model, num_peaks):
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

        mag_model: the type of magnetism the peak corresponds to

        num_peaks: the number of peaks in the system

    Return:

        intesity: the neutron intensity

    """

    intensity = 0

    for i in range(num_peaks):

        if(mag_model[i] == 'FM'):

            intensity += -A*chi_imag_FM_func(coords, chi_Q_inv[i], gamma[i], c[i], Q[i])

        elif(mag_model[i] == 'AFM'):

            intensity += -A*chi_imag_AFM_func(coords, chi_Q_inv[i], gamma[i], c[i], Q[i])

        else:

            print("Wrong inputs in the magnetic model list, please select either FM or AFM")
            exit()

    return intensity


def wrapper_intensity_func(coords, mag_model, num_peaks, *args):
    """
    A wrapper function to convert the input for the intensity function 
    data from a list of all params in args to a lists for each parameter

    Args:

        coords: an 2d array containing the q space coords and the frequency 
        coords order is [q,w]

        mag_model: the type of magnetism the peak corresponds to

        num_peaks: the number of peaks in the system

        args: system parameters in single list format

    Return:

        intesity: the neutron intensity

    """
    #extract variable parameters
    A = args[0][0]
    chi_Q_inv, gamma, c, Q = list(args[0][1:num_peaks + 1]), list(args[0][num_peaks + 1: 2*num_peaks + 1]), list(args[0][2*num_peaks + 1: 3*num_peaks + 1]), list(args[0][3*num_peaks + 1: 4*num_peaks + 1])

    #calculate intensity
    intensity = intensity_func(coords, A, chi_Q_inv, gamma, c, Q, mag_model, num_peaks)

    return intensity


#main script
if __name__ == "__main__":

    #Data params
    num_peaks = 2
    #peak types
    mag_model = ['FM', 'AFM']

    #initial params
    A = [-1.0]
    chi_Q_inv = [1.0, 1.0]
    gamma = [1.0, 1.0]
    c = [1000.0, 1000.0]
    Q = [0.25, 0.75]

    #create parameters list
    params = A + chi_Q_inv + gamma + c + Q

    q = np.linspace(0, 1, 1000)
    freq = np.full((1000),1.0)

    coords = np.array([q,freq], dtype=object)

    coords = coords.T

    intensity = wrapper_intensity_func(coords, mag_model, num_peaks, params)

    plt.plot(q, intensity)

    plt.show()

    #fit curve
    #popt, pcov = optimize.curve_fit(lambda x, *params_init: wrapper_intensity_func(coords, mag_model, num_peaks, params_init), coords_in, data, p0=params_init)









    








