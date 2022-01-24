from contextvars import copy_context
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import scipy.cluster.vq as scv
import math 

def Colour_To_Array(arr,cmap):  
    """
    Convert a colourmapped image to an array of scalar values corresponding the the values 
    at that point

    Args:

        arr: the RGB(A) converted array

        cmap: the colourmap being used

    Returns:

        values: conversion to scalar values
    """  
    # convert colourmap to a gradient of scalar values
    gradient = cmap(np.linspace(0.0,1.0,10000))

    # Reshape input array to something like (N*N, 3/4), where N is number of pixels in input and 
    # the third index is a fourtuple of RGB(A) values
    arr2 = arr.reshape((arr.shape[0]*arr.shape[1],arr.shape[2]))

    #change gradient values depending on whether RGBA or RGB
    gradient = gradient[:,0:arr.shape[2]]

    # Use vector quantization to shift the values in arr2 to the nearest point in
    # the code book (gradient).
    code, _ = scv.vq(arr2,gradient)

    # code is an array of length arr2 (N, N), holding the code book index for
    # each observation. (arr2 are the "observations".)
    # Scale the values so they are from 0 to 1.
    values = code.astype('float')/gradient.shape[0]

    #replace values with NaN
    values[values==code.astype('float')[0]/gradient.shape[0]] = np.nan

    # Reshape values back to (N,N)
    values = values.reshape(arr.shape[0],arr.shape[1])
    values = values[::-1]
    return values


def Convert_To_Scaled(values, e_min, e_max, e_interval, q_min, q_max, q_interval):
    """
    Convert array of values to scaled data points for fitting 

    Args:

        values: array containing data in 2D array format

        e_min: the min energy range value

        e_max: the max energy range value

        e_interval: the step at which to sample the image in energy space

        q_min: the min q range value

        q_max: the max q range value

        q_interval: the step at which to sample the image in q space

    Returns:

        coord_data: data in coordinate form for fitting
    """

    height, width = values.shape

    e_scale = np.linspace(e_min, e_max, height)
    e_step = int(height*e_interval/(e_max - e_min))

    q_scale = np.linspace(q_min, q_max, width)
    q_step = int(width*q_interval/(q_max - q_min))

    coord_arrays = []

    for i in range(math.floor((height)/e_step) + 1):
        for j in range(math.floor((width)/q_step) + 1):

            h_index = i*e_step
            w_index = j*q_step

            temp_val = [q_scale[w_index], e_scale[h_index], values[h_index, w_index]]

            if not(math.isnan(values[h_index, w_index])):
                coord_arrays += (temp_val,)

    coord_data = np.array(coord_arrays)

    return coord_data


    

arr = plt.imread('data.png')
values = Colour_To_Array(arr,cm.turbo)

out_data = Convert_To_Scaled(values, 2, 16, 0.2, 0.3, 1.2, 0.1)

np.savetxt("out_data.csv", out_data, delimiter=",")
# Proof that it works:
plt.imshow(values,interpolation='bilinear', cmap=cm.turbo,
           origin='lower', extent=[-3,3,-3,3])
plt.show()