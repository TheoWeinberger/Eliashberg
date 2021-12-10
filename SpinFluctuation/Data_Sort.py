#--------------------------------------------------------------------------
# Data_Sort.py - Python package to sort NR data from ISIS for YFe2Ge2
#
# version: 1
# date: 09/12/2021
# author: Theo Weinberger
#--------------------------------------------------------------------------

from os import error
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def Clean_Data(filename):
    """
    Function to clean the raw data from the ISIS data files
    and reduce it down to the minimal dataset

    Args: 

        filename: string containing the name of the file 
        being cleaned

    Return:

        df_out: a dataframe containing the cleaned data

    """
    #read in raw data
    df = pd.read_csv(filename, sep = " ", engine="python", header=None)

    #list to store cleaned data
    df_whole_list = []

    #iterate though raw data and remove any nan values
    #this produces the cleaned data
    for j in range(df.shape[0]):

        df_row_list = []

        for i in range(df.shape[1]):

            if pd.isnull(df[i].iloc[j]) == False:

                df_row_list += [df[i].iloc[j]]

        df_whole_list += [df_row_list]

    #create the new cleaned dataframe
    df_clean = pd.DataFrame.from_records(df_whole_list)

    #calculate the intensity as counts/second
    df_clean["Intensity"] =  df_clean[8]/df_clean[7]
    #and the error as poissonian error
    df_clean["Error"] = df_clean["Intensity"]/10.0

    #name column names for reduced output
    column_names = ["Q", "Energy", "Intensity", "Error"]

    #create output datafrane
    df_out = pd.DataFrame(columns=column_names)

    ##assign new columns with data rounded
    df_out["Q"] = df_clean[1].round(2)
    df_out["Energy"] = df_clean[4].round(2)
    df_out["Intensity"] = df_clean["Intensity"]
    df_out["Error"] = df_clean["Error"]

    return df_out


def Normalise_Q(df_q, q_norm):
    """
    Function to normalise data inputs that consist of Q varying data

    Args:

        df_q: dataframe containing the q dependent data

        q_norm: value of q against which the data should be normalised

    Returns:

        df_q_norm: renormalised dataframe
    """
    #initialise normed data
    df_q_norm = df_q

    #get normalisation values
    norm_value = df_q_norm[df_q_norm["Q"] == q_norm]["Intensity"]
    error_norm = df_q_norm[df_q_norm["Q"] == q_norm]["Error"]

    #calculate renormalised values
    df_q_norm["Intensity"] -=  float(norm_value)
    df_q_norm["Error"] =  np.sqrt(df_q_norm["Error"]**2 + float(error_norm)**2)

    return df_q_norm

def Normalise_E(df_e, df_e_base):
    """
    Function to normalise data inputs that consist of energy varying data

    Args:

        df_e: dataframe containing the e dependent data

        df_e_base: dataframe containing the e dependent data measured
        at the background reading point

    Returns:

        df_e_norm: renormalised dataframe
    """
    #initialise data
    df_e_norm =  df_e

    #calculate normed data    
    df_e_norm["Intensity"] =  df_e["Intensity"] - df_e_base["Intensity"]
    df_e_norm["Error"] =  np.sqrt(df_e["Error"]**2 + df_e_base["Error"]**2)

    return df_e_norm

def Clean_Energy_Data(energy_data, background_data):
    """
    Function to clean the raw data from the ISIS data files
    and reduce it down to the minimal dataset for the case 
    of energy varying data

    Args: 

        energy_data: string containing the name of the file 
        being cleaned

        background_data: string containing the name of the file 
        being used as background reference

    Return:

        df_out: a dataframe containing the cleaned and normalised
        data

    """
    #load in files and clean
    df_e = Clean_Data(energy_data)
    df_e_base = Clean_Data(background_data)

    #normalise
    df_out = Normalise_E(df_e, df_e_base)

    return df_out


def Clean_Q_Data(q_data, q_norm):
    """
    Function to clean the raw data from the ISIS data files
    and reduce it down to the minimal dataset for the case 
    of q varying data

    Args: 

        q_data: string containing the name of the file 
        being cleaned

        q_norm: the q value being used as the background 
        reading

    Return:

        df_out: a dataframe containing the cleaned and normalised
        data

    """
    #load in files and clean
    df_q = Clean_Data(q_data)
    #normalise
    df_out = Normalise_Q(df_q, q_norm)

    return df_out



#main script
if __name__ == "__main__":

    #set these variables depending on what you
    #want to normalise

    #name of input file
    filename1 = "YFe2Ge2QE35.dat"
    filename2 = "YFe2Ge2QE35.dat"

    #name of output file
    fileout = "YFe2Ge2QE35_Clean.dat"

    #whether it's q or e data being cleaned    
    clean_type = "q"

    #q value being used as background
    q_norm = 0.3

    if clean_type == "q":

        df_out = Clean_Q_Data(filename1, q_norm)

    elif clean_type == "e":

        df_out = Clean_Energy_Data(filename1, filename2)

    else:

        print("Invalid value for clean_type. Please choose either q or e.")


    df_out.to_csv(fileout)






