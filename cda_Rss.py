"""
 *========================================================================================
 *
 *       File name:  cda_Rss.py
 *
 *    Description:  Takes Q vs T data collected during a cooldown and  plots Rs* vs T
 *
 *        Created:  27/11/2021
 *       Compiler:  Python 3.8
 *
 *         Author:  Ruth Gregory
 *          Email:  ruthannrgregory@gmail.com
 *   Organization:  TRIUMF
 *
 * =======================================================================================
 """
print("Running the script")
print()

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
import numpy as np
import pandas as pd

"""---------------------------------------------------------------------------------------
Get all the input parameters ready
---------------------------------------------------------------------------------------"""
#These are the beta values to use in the fit function corrections when converting Rs* to Rs
#The six entries in the lists correspond to the six different coaxial cavities
beta0_all = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

beta1_all = [1.43265730367638, 1.4731736554568,
    1.46267321592333, 1.46116617646345,
    1.46289120717899, 1.46313972701761]

beta2_all = [1.77803816661436, 1.8706480191515,
    1.85703995599432, 1.85657359333326,
    1.82469356803503,1.86243947537751]

beta3_all = [2.06107129508082, 2.21325662259346,
    2.19252463157984, 2.1964531414152,
    2.20609678777995, 2.20709554483536]

#Possible geometric factor and frequency values for the six cavities
G_vals = [37.47, 113.7, 60.39, 120.77, 181.08, 241.24]
frequencies = [217, 647, 389, 778, 1166, 1555]

#Get file path and name from user
file_path = input("Enter the full path and name of your input data file: ")
file = open(file_path, "r")
"""
#Get information about the input data table from user
Qcol_str = input("Enter the column number of your Quality Factor data: ")
Tcol_str = input("Enter the column number of your temperature data: ")
Eacccol_str = input("Enter the column number of your accelerating field data: ")
field_vals_str = input("Enter the approximate target accelerating field values separated by a space: ")
skiplines_str = input("Optional input: Enter the number of lines to skip at the beginning of the table, or to not skip any lines, press enter: ")

msg = (
    '{sep}'
    'Enter a number to choose your coaxial cavity: {sep}'
    '{sep}'
    'Enter 1 for QWR 217 MHz with G = 37.47 Ohms {sep}'
    'Enter 2 for QWR 648 MHz with G = 113.7 Ohms {sep}'
    'Enter 3 for HWR 389 MHz with G = 60.39 Ohms {sep}'
    'Enter 4 for HWR 778 MHz with G = 120.77 Ohms {sep}'
    'Enter 5 for HWR 1166 MHz with G = 181.08 Ohms {sep}'
    'Enter 6 for HWR 1555 MHz with G = 241.24 Ohms'
    '{sep}').format(sep='\n')
print(msg)

in_cavity = input("Enter number: ") #Value for gemetric factor, freqeuncy, and determining which betas to use
"""

skiplines = 0 #The default number of lines to skip at the beginning of the table is zero
"""
#Get the list of accelerating field values
FieldValues = []
field_vals_list = field_vals_str.split()
for val in field_vals_list:
    FieldValues.append(float(val))

#Convert the user input values to ints and get parameters
if skiplines_str != "":
    skiplines = int(skiplines_str)

Qcol = int(Qcol_str)
Tcol = int(Tcol_str)
Eacccol = int(Eacccol_str)
cavity = int(in_cavity)
SWR = 1.080124
""""""
#Subtract 1 because Python indexing starts at 0
Qcol -= 1
Tcol -= 1
Eacccol -= 1
cavity -= 1
"""
SWR = 1.20534
FieldValues = [10,20,30,40,50,60,70,80,90]
Qcol = 10
Tcol = 21
Eacccol = 8
cavity = 1
#Get frequency and G factor values
freq = frequencies[cavity]
G = G_vals[cavity]

#Create empty lists to hold the data
Qdata = []
Tdata = []
Eaccdata = []
Rss = []

Bp_err = []
Rss_err = []
weights = []

"""---------------------------------------------------------------------------------------
Read in the data from the input file
---------------------------------------------------------------------------------------"""
#Read the input file line by line and extract data
line_cnt = 1 #line count, it will start with line 1
for line in file:

    if 'New calibration' in line:
        file.readline()
        file.readline()

    if line_cnt <= skiplines: #skip the indicated number of lines at the beginning
        line_cnt += 1
        continue

    if line_cnt == 1:
        line = line.replace(' ','')
        columns = line.split()
        col_cnt = 0

        for label in columns:
            if 'TEMPK3' in label:
                Tcol = col_cnt
            col_cnt += 1

    line = line.strip()
    columns = line.split()

    if columns == []: #skip any columns that are completely empty
        line_cnt += 1
        continue

    try: #add data to the lists
        if float(columns[Qcol]) > 0:
            Qdata.append(float(columns[Qcol]))
        else:
            continue
    except IndexError:
        print("Warning: index error on line", line_cnt, "in the Q column. Skipping line", line_cnt)
        line_cnt += 1
        continue
    except ValueError:
        print("Value on line", line_cnt, "in the Q column is not a number. Skipping line", line_cnt)
        line_cnt += 1
        continue

    try:
        Eaccdata.append(float(columns[Eacccol]))
    except IndexError:
        print("Warning: index error on line", line_cnt, "in the Eacc column. Skipping line", line_cnt)
        Qdata.pop()
        line_cnt += 1
        continue
    except ValueError:
        print("Value on line", line_cnt, "in the Eacc column is not a number. Skipping line", line_cnt)
        Qdata.pop()
        line_cnt += 1
        continue

    try: #Skip an entire line if any of the columns has an error so that the lists are of the same length
        if float(columns[Tcol]) > 0:
            Tdata.append(float(columns[Tcol]))
        else:
            Tdata.append(np.nan)
            line_cnt += 1
            continue
    except IndexError:
        Tdata.append(np.nan)
        line_cnt += 1
        continue
    except ValueError:
        print("Value on line", line_cnt, "in the Temperature column is not a number. Skipping line", line_cnt)
        Qdata.pop()
        Eaccdata.pop()
        line_cnt += 1
        continue

    line_cnt += 1

file.close()
"""---------------------------------------------------------------------------------------
For each accelerating field ramp up, calculate the feild corrected Rs values
---------------------------------------------------------------------------------------"""
Tdata_pd = pd.Series(Tdata)
Tdata = Tdata_pd.interpolate()
#print(Tdata.to_string())

#Get the non-corrected Rs* data
for i in range(0,len(Qdata)):
    Rss.append((G/Qdata[i])*10**9)
    Bp_err.append(Eaccdata[i]*(SWR-1)/4)
    Rss_err.append(Rss[i]*(SWR-1)/2)
    weights.append(1 - (Rss_err[i]/Rss[i]))

"""---------------------------------------------------------------------------------------
Separate the data by different field amplitude values
---------------------------------------------------------------------------------------"""

Inv_Tdata_sep = [] #list of lists for inverse temperature data separated by field amplitude
Tdata_sep = [] #list of listst for temperature data separated by field amplitude
Rss_sep = [] #list of lists for surface resistance data in nano-ohms separated by field amplitude
Rss_err_sep = []
Rss_sep_log = [] #list of lists for the logarithm of the surface resistance data in nano-ohms separated by field amplitude
Rss_sep_ln = [] #list of lists for the natural logarithm of the surface resistance data in nano-ohms separated by field amplitude
weights_sep = []

for value in FieldValues:
    Inv_Tdata_sep.append([])
    Tdata_sep.append([])
    Rss_sep.append([])
    Rss_err_sep.append([])
    Rss_sep_log.append([])
    Rss_sep_ln.append([])
    weights_sep.append([])

#Make two lists containing lists of the inverse temperature data and the
#corrected surface resistance data. Each sub-list corresponds to a different
#field amplitude

for i in range(0,len(Eaccdata)):

    for j in range(0,len(FieldValues)):

        if (Eaccdata[i] < FieldValues[j]+0.5) and (Eaccdata[i] > FieldValues[j]-0.5):
            Inv_Tdata_sep[j].append(1/Tdata[i])
            Tdata_sep[j].append(Tdata[i])
            Rss_sep[j].append(Rss[i])
            Rss_err_sep[j].append(Rss_err[i])
            Rss_sep_log[j].append(np.log10(Rss[i]))
            Rss_sep_ln[j].append(np.log(Rss[i]))
            weights_sep[j].append(weights[i])

"""---------------------------------------------------------------------------------------
Plot Rs* vs T
---------------------------------------------------------------------------------------"""
fig1, ax1 = plt.subplots(nrows=1, ncols=1)

colors = ['b','orange', 'g', 'r', 'c', 'm', 'y', 'salmon', 'brown', 'lawngreen' , '0.4', '0.8' ]
shapes = ['^', 's', 'P', '*', '+', 'd', 'x']

legend_entries = []
temps_legend = []

#Temperature values for the fit plotting
T_vals = np.linspace(1.9, 4.5, 50)
Inv_T_vals = np.linspace(1/1.9, 1/4.5, 50)

for i in range(0,len(FieldValues)):

    name = str(FieldValues[i]) + " mT"
    legend_entries.append(name)

if len(legend_entries) > (len(colors)+len(shapes)):
    print("Warning: This program wasn't expecting more than 19 different field amplitudes. Some field amplitude data sets will be indistiguishable on the plot, as they will be plotted as black dots.")


for i in range(0,len(legend_entries)):

    if i < len(colors):
        #ax1.errorbar(Inv_Tdata_sep[i],Rs_sep[i], yerr=Rs_err_sep[i], marker='o', linestyle='none', markersize=4, color=colors[i], label=legend_entries[i])
        ax1.plot(Inv_Tdata_sep[i],Rss_sep[i], marker='o', linestyle='none', markersize=4, color=colors[i], label=legend_entries[i])

    elif i < (len(colors)+len(shapes)):
        ax1.plot(Inv_Tdata_sep[i],Rs_sep[i], marker=shapes[i-len(colors)], linestyle='none', markersize=4, color='black', label=legend_entries[i])
        ax1.plot(Inv_Tdata_sep[i], result.best_fit, marker='None', linestyle='--',color='black')
        ax2.plot(Inv_Tdata_sep[i],((result.residual))*100/Rs_sep[i], marker=shapes[i-len(colors)], markersize=3, color='black', label=legend_entries[i])
    else:
        ax1.plot(Inv_Tdata_sep[i],Rs_sep[i], marker='o', linestyle='none', markersize=4, color='black', label=legend_entries[i])
        ax1.plot(Inv_Tdata_sep[i], result.best_fit, marker='None', linestyle='--',color='black')
        ax2.plot(Inv_Tdata_sep[i],((result.residual))*100/Rs_sep[i], marker='o', markersize=3, color='black', label=legend_entries[i])


#format the x axis
x_formatter = FixedFormatter([r'10$^{-1}$', r'9$^{-1}$', r'8$^{-1}$', r'7$^{-1}$', r'6$^{-1}$', r'5$^{-1}$',
 r'4.5$^{-1}$', r'4.0$^{-1}$', r'3.5$^{-1}$', r'3.0$^{-1}$', r'2.5$^{-1}$', r'2.2$^{-1}$', r'2.0$^{-1}$', r'1.8$^{-1}$', r'1.7$^{-1}$'])
x_locator = FixedLocator([0.1, 1/9, 0.125, 1/7, 1/6, 0.2, 1/4.5, 0.25, 1/3.5, 1/3, 0.4, 1/2.2, 0.5, 1/1.8, 1/1.7])

#ax1.set_yscale('log')
ax1.set_xlabel(r'Inverse Temperature [K$^{-1}$]',fontsize=14)
ax1.set_ylabel(r'R$_s*[n\Omega]$',fontsize=14)
ax1.xaxis.set_major_formatter(x_formatter)
ax1.xaxis.set_major_locator(x_locator)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.vlines(x=0.45977, ymin=0, ymax=200, color='gray', linestyle='dashed', label = 'lambda point')
ax1.legend(title='Field Amplitude',fontsize=14)
ax1.set_title("Surface Resistance vs Inverse Temperature",fontsize=18)
ax1.grid(True)


plt.show()
