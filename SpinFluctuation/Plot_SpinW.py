#--------------------------------------------------------------------------
# Plot_SpinW.py - Python script to plot spinW data against ILL YFe2Ge2 data
#
# version: 1
# date: 08/02/2022
# author: Theo Weinberger
#--------------------------------------------------------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

cp = sns.color_palette("bright")

plt.style.use(['science', 'grid'])


def find_nearest(array, value):
    """
    Function to find index of item in an array closest
    to the desired value

    Args:

        array: array to search

        value: value to find

    Returns:

        idx: array of indeces of points in array closest to this value
    """
    array = np.asarray(array)
    idx = np.nanargmin((np.abs(array - value)))

    return idx


#First import YFG SpinW simulation data
spinw_data = np.genfromtxt('data/YFe2Ge2_SpinW.dat', delimiter=',')

spinw_data = 1.56*spinw_data

#define axes, define the beginning and end points for the data here
e_range = np.linspace(0, 60, spinw_data.shape[0])
q_range = np.linspace(-0.3, 1.0, spinw_data.shape[1])

#plot data
qq, ee = np.meshgrid(q_range, e_range)
plt.set_cmap('jet')
plt.contourf(qq, ee, spinw_data)
plt.ylabel("Energy/meV")
plt.xlabel("hh3/r.l.u")
plt.ylim(0,16)
plt.xlim(-0.2,0.8)
plt.savefig('EHHK_SW_Global.pdf')
plt.show()


#cut along q axis at FM peak
q_cut = find_nearest(q_range, 0)
plot_data = spinw_data[:,q_cut]

#comparison data
df_plot = pd.read_csv("data/YFe2Ge2EQ0_clean_m1.dat", sep = ",", engine="python")

#protect against divide by zero
#errors_plot = np.divide(errors_plot, abs(data), out=np.zeros_like(errors_plot), where=data>=0.1)
#errors_plot[errors_plot == 0] = np.partition(np.unique(errors_plot),1)[1]
plt.plot(e_range, plot_data, label="Fit", color='black', linestyle='--')
plt.errorbar(df_plot["Energy"], df_plot["Intensity"], df_plot["Error"], label="Data", ls='none', marker='+', color='black', capsize=2)
plt.legend(prop={'size': 6})
plt.xlabel("Energy/meV")
plt.ylabel("Normalised Counts")

plt.xlim(0,14)
plt.savefig('AFM_Energy_SW_Global.pdf')
plt.show()


#cut along q axis at AFM peak
q_cut = find_nearest(q_range, 0.5)
plot_data = spinw_data[:,q_cut]

#comparison data
df_plot = pd.read_csv("data/YFe2Ge2EQ0_5_clean_m1.dat", sep = ",", engine="python")

#protect against divide by zero
#errors_plot = np.divide(errors_plot, abs(data), out=np.zeros_like(errors_plot), where=data>=0.1)
#errors_plot[errors_plot == 0] = np.partition(np.unique(errors_plot),1)[1]
plt.plot(e_range, plot_data, label="Fit", color='black', linestyle='--')
plt.errorbar(df_plot["Energy"], df_plot["Intensity"], df_plot["Error"], label="Data", ls='none', marker='+', color='black', capsize=2)
plt.legend(prop={'size': 6})
plt.xlabel("Energy/meV")
plt.ylabel("Normalised Counts")

plt.xlim(0,14)
plt.savefig('FM_Energy_SW_Global.pdf')
plt.show()


#cut along energy axis at 6 meV
e_cut = find_nearest(e_range, 6)
plot_data = spinw_data[e_cut,:]

#comparison data
df_plot = pd.read_csv("data/YFe2Ge2QE6_clean_m1.dat", sep = ",", engine="python")

#protect against divide by zero
#errors_plot = np.divide(errors_plot, abs(data), out=np.zeros_like(errors_plot), where=data>=0.1)
#errors_plot[errors_plot == 0] = np.partition(np.unique(errors_plot),1)[1]
plt.plot(q_range, plot_data, label="Fit", color='black', linestyle='--')
plt.errorbar(df_plot["Q1"], df_plot["Intensity"], df_plot["Error"], label="Data", ls='none', marker='+', color='black', capsize=2)
plt.legend(prop={'size': 6})
plt.xlabel("q/r.l.u")
plt.xlim(-0.2,0.8)
plt.ylabel("Normalised Counts")
plt.savefig('E6_QHH3_SW_Global.pdf')
plt.show()


#cut along energy axis at 3.5 meV
e_cut = find_nearest(e_range, 3.5)
plot_data = spinw_data[e_cut,:]

#comparison data
df_plot = pd.read_csv("data/YFe2Ge2QE3_5_clean_m1.dat", sep = ",", engine="python")

#protect against divide by zero
#errors_plot = np.divide(errors_plot, abs(data), out=np.zeros_like(errors_plot), where=data>=0.1)
#errors_plot[errors_plot == 0] = np.partition(np.unique(errors_plot),1)[1]
plt.plot(q_range, plot_data, label="Fit", color='black', linestyle='--')
plt.errorbar(df_plot["Q1"], df_plot["Intensity"], df_plot["Error"], label="Data", ls='none', marker='+', color='black', capsize=2)
plt.legend(prop={'size': 6})
plt.xlabel("q/r.l.u")
plt.xlim(-0.2,0.8)
plt.ylabel("Normalised Counts")
plt.savefig('E3_5_QHH3_SW_Global.pdf')
plt.show()

#cut along energy axis at 3.5 meV
e_cut = find_nearest(e_range, 2)
plot_data = spinw_data[e_cut,:]

#comparison data
df_plot = pd.read_csv("data/YFe2Ge2QE2_clean_m1.dat", sep = ",", engine="python")

#protect against divide by zero
#errors_plot = np.divide(errors_plot, abs(data), out=np.zeros_like(errors_plot), where=data>=0.1)
#errors_plot[errors_plot == 0] = np.partition(np.unique(errors_plot),1)[1]
plt.plot(q_range, plot_data, label="Fit", color='black', linestyle='--')
plt.errorbar(df_plot["Q1"], df_plot["Intensity"], df_plot["Error"], label="Data", ls='none', marker='+', color='black', capsize=2)
plt.legend(prop={'size': 6})
plt.xlabel("q/r.l.u")
plt.xlim(-0.2,0.8)
plt.ylabel("Normalised Counts")
plt.savefig('E2_QHH3_SW_Global.pdf')
plt.show()


fig, ax = plt.subplots()

plt.plot(q_range, plot_data, label="2 meV",  linestyle='--', color = cp[0])
plt.errorbar(df_plot["Q1"], df_plot["Intensity"], df_plot["Error"], label="2 meV", ls='none', marker='+',  capsize=2, color = cp[0])

df_plot = pd.read_csv("data/YFe2Ge2QE3_5_clean_m1.dat", sep = ",", engine="python")
#cut along energy axis at 3.5 meV
e_cut = find_nearest(e_range, 3.5)
plot_data = spinw_data[e_cut,:]
#protect against divide by zero
#errors_plot = np.divide(errors_plot, abs(data), out=np.zeros_like(errors_plot), where=data>=0.1)
#errors_plot[errors_plot == 0] = np.partition(np.unique(errors_plot),1)[1]
plt.plot(q_range, plot_data, label="3.5 meV",  linestyle='--', color = cp[1])
plt.errorbar(df_plot["Q1"], df_plot["Intensity"], df_plot["Error"], label="3.5 meV", ls='none', marker='+',  capsize=2, color = cp[1])

#comparison data
df_plot = pd.read_csv("data/YFe2Ge2QE6_clean_m1.dat", sep = ",", engine="python")
#cut along energy axis at 3.5 meV
e_cut = find_nearest(e_range, 6)
plot_data = spinw_data[e_cut,:]
#protect against divide by zero
#errors_plot = np.divide(errors_plot, abs(data), out=np.zeros_like(errors_plot), where=data>=0.1)
#errors_plot[errors_plot == 0] = np.partition(np.unique(errors_plot),1)[1]
plt.plot(q_range, plot_data, label="6 meV", linestyle='--', color = cp[2])
plt.errorbar(df_plot["Q1"], df_plot["Intensity"], df_plot["Error"], label="6 meV", ls='none', marker='+',  capsize=2, color = cp[2])

plt.legend(prop={'size': 6}, ncol=2)
plt.xlabel("q/r.l.u")
plt.xlim(-0.2,0.8)
plt.ylabel("Normalised Counts")
h, l = ax.get_legend_handles_labels()
ph = [plt.plot([],marker="", ls="")[0]]*2
handles = ph[:1] + [h[::2][1]] + [h[1::2][0]] + [h[::2][0]] + ph[1:] +  [h[1::2][2]] + [h[::2][2]] + [h[1::2][1]] 


labels = ["Fit"] + [l[::2][1]] + [l[1::2][0]] + [l[::2][0]] + ["Data"] + [l[1::2][2]] + [l[::2][2]] + [l[1::2][1]] 
leg = plt.legend(handles, labels, ncol=2, prop={'size': 6})

for vpack in leg._legend_handle_box.get_children():
    for hpack in vpack.get_children()[:1]:
        hpack.get_children()[0].set_width(0)
plt.savefig('QHH3_SW_Global.pdf')
plt.show()


