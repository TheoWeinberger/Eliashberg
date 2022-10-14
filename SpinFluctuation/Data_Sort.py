#--------------------------------------------------------------------------
# Data_Sort.py - Python package to sort NR data from ISIS for YFe2Ge2
#
# version: 1
# date: 09/12/2021
# author: Theo Weinberger
#--------------------------------------------------------------------------

from os import error
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams

def Clean_Data(filename, norm):
    """
    Function to clean the raw data from the ISIS data files
    and reduce it down to the minimal dataset

    Args: 

        filename: string containing the name of the file 
        being cleaned

        norm: The normalisation reference

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

    if norm == "time":
    
        #calculate the intensity as counts/second
        df_clean["Intensity"] =  df_clean[8]/df_clean[7]
        #and the error as poissonian error
        df_clean["Error"] = np.sqrt(df_clean[8])/df_clean[7]
    
    elif norm == "m1":
        
        #calculate the intensity as 10,000*counts/m1
        df_clean["Intensity"] =  10000*df_clean[8]/df_clean[5]
        #and the error as poissonian error
        df_clean["Error"] = 10000*np.sqrt(df_clean[8])/df_clean[5]

    elif norm == "m2":
        
        #calculate the intensity as counts/m2
        df_clean["Intensity"] =  df_clean[8]/df_clean[6]
        #and the error as poissonian error
        df_clean["Error"] = np.sqrt(df_clean[8])/df_clean[6]

    else:

        print("Invalid value for normalisation type, choose time, m1 or m2")
        exit(1)
    
    

    df_clean["Temp"] = df_clean[9]

    #name column names for reduced output
    column_names = ["Q1", "Q2", "Q3", "Energy", "Intensity", "Error", "Temp"]

    #create output datafrane
    df_out = pd.DataFrame(columns=column_names)

    ##assign new columns with data rounded
    df_out["Q1"] = df_clean[1].round(2)
    df_out["Q2"] = df_clean[2].round(2)
    df_out["Q3"] = df_clean[3].round(2)
    df_out["Energy"] = df_clean[4].round(2)
    df_out["Intensity"] = df_clean["Intensity"]
    df_out["Error"] = df_clean["Error"]
    df_out["Temp"] = df_clean["Temp"]


    #remove 1/q^2 dependence
    #df_out["Intensity"] = df_out["Intensity"]*np.sqrt(df_out["Q1"]**2 + df_out["Q2"]**2 + df_out["Q3"]**2)

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
    norm_value = df_q_norm[df_q_norm["Q3"] == q_norm]["Intensity"]
    error_norm = df_q_norm[df_q_norm["Q3"] == q_norm]["Error"]
    temp_norm = df_q_norm[df_q_norm["Q3"] == q_norm]["Temp"]
    energy_norm = df_q_norm[df_q_norm["Q3"] == q_norm]["Energy"]

    boltz = 8.617e-5

    #data must be scaled for this
    bose_norm = 1.0/(1.0 - np.exp(-(float(energy_norm)*1e-3)/(boltz*float(temp_norm))))

    #calculate renormalised values
    df_q_norm["Intensity"] = df_q_norm["Intensity"]/Bose(df_q_norm) - float(norm_value)/bose_norm
    df_q_norm["Error"] =  np.sqrt((df_q_norm["Error"]/Bose(df_q_norm))**2 + float(error_norm/bose_norm)**2)

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
    df_e_norm["Error"] =  np.sqrt((df_e["Error"]/Bose(df_e))**2 + (df_e_base["Error"]/Bose(df_e_norm))**2)

    return df_e_norm

def Clean_Energy_Data(energy_data, background_data, norm, threshold):
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
    df_e = Clean_Data(energy_data, norm)
    df_e_base = Clean_Data(background_data, norm)


    df_e = df_e.sort_values(by=["Energy"])
    df_e_base = df_e_base.sort_values(by=["Energy"])

    #apply bose corrections
    df_e["Intensity"] = df_e["Intensity"]/Bose(df_e)
    df_e_base["Intensity"] = df_e_base["Intensity"]/Bose(df_e_base)

    #apply q scaling to normalisation data
    scaling = (df_e_base["Q1"]**2 + df_e_base["Q2"]**2 + df_e_base["Q3"]**2)/(df_e["Q1"]**2 + df_e["Q2"]**2 + df_e["Q3"]**2)

    df_e_base["Intensity"] = df_e_base["Intensity"]*scaling

    df_e["Error"] = df_e["Error"]/Bose(df_e)
    df_e_base["Error"] = df_e_base["Error"]/Bose(df_e_base)

    #crop data to 2meV and above
    df_e = df_e.drop(df_e[df_e["Energy"] < threshold].index)
    df_e_base = df_e_base.drop(df_e_base[df_e_base["Energy"] < threshold].index)

    #linear fitting to background data 
    energy = df_e_base["Energy"]
    counts = df_e_base["Intensity"]
    slope, intercept, r_value, p_value, std_err = stats.linregress(energy, counts)

    df_e_base_fit = df_e.copy()

    fit_values = df_e["Energy"]*slope + intercept

    df_e_base_fit["Intensity"] = fit_values

    #plot raw data and subsequent fittings

    plt.errorbar(df_e["Energy"], df_e["Intensity"], df_e["Error"], ls='none', marker='x', color='black', capsize=2)
    plt.errorbar(df_e_base["Energy"], df_e_base["Intensity"], df_e_base["Error"], ls='none', marker='+', color='black', capsize=2)
    plt.plot(df_e_base_fit["Energy"], df_e_base_fit["Intensity"], color='black', linestyle='--')
    plt.xlabel("Energy/meV")
    plt.ylabel("Normalised Counts")
    plt.grid(linestyle='--', linewidth='0.5', color='gray')
    plt.show()



    #match data points in df_e to df_e_base

    #normalise
    df_out = Normalise_E(df_e, df_e_base_fit)



    return df_out


def Clean_Q_Data(q_data, q_norm, norm):
    """
    Function to clean the raw data from the ISIS data files
    and reduce it down to the minimal dataset for the case 
    of q varying data

    Args: 

        q_data: string containing the name of the file 
        being cleaned

        q_norm: the q value being used as the background 
        reading

        norm: The normalisation reference

    Return:

        df_out: a dataframe containing the cleaned and normalised
        data

    """
    #load in files and clean
    df_q = Clean_Data(q_data, norm)
    #normalise
    df_out = Normalise_Q(df_q, q_norm)

    df_out.sort_values(by=["Q3"])

    return df_out


def Bose(data):
    """
    Calculate the bose function for emission for 
    a data frame

    Args:

        data: a dataframe containing scattering data

    Return

        bose: bose factors
    """

    boltz = 8.617e-5

    #data must be scaled for this
    bose = 1.0/(1.0 - np.exp(-(data["Energy"]*1e-3)/(boltz*data["Temp"])))

    return bose



#main script
if __name__ == "__main__":

    #set these variables depending on what you
    #want to normalise

    #name of input file
    filename1 = "data/YFe2Ge2L1E_6.dat"
    filename2 = "data/YFe2Ge2E_back.dat"

    #name of output file
    fileout = "data/YFe2Ge2L1E_6_clean_m1.dat"

    #whether it's q or e data being cleaned    
    clean_type = "q"

    #Normalisation type, time, m1 or m2
    norm = "m1"

    #q value being used as background
    q_norm = 0.5

    if clean_type == "q":

        df_out = Clean_Q_Data(filename1, q_norm, norm)

        plt.errorbar(df_out["Q3"], df_out["Intensity"], df_out["Error"], ls='none', marker='+', color='black', capsize=2)
        plt.xlabel("r.l.u.")
        plt.ylabel("Normalised Counts")
        plt.grid(linestyle='--', linewidth='0.5', color='gray')
        plt.show()

    elif clean_type == "e":

        df_out = Clean_Energy_Data(filename1, filename2, norm, 1)

        plt.errorbar(df_out["Energy"], df_out["Intensity"], df_out["Error"], ls='none', marker='+', color='black', capsize=2)
        plt.xlabel("Energy/meV")
        plt.ylabel("Normalised Counts")
        plt.grid(linestyle='--', linewidth='0.5', color='gray')
        plt.show()

    else:

        print("Invalid value for clean_type. Please choose either q or e.")


    df_out.to_csv(fileout)






