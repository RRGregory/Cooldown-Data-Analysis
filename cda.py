"""
 *========================================================================================
 *
 *       File name:  cda.py
 *
 *    Description:  Takes Q vs T data collected during a cooldown and extracts
 *                  Rbcs and Rres as a function of field amplitude
 *
 *        Created:  19/05/2021
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
from scipy.optimize import curve_fit
from lmfit import Model

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

fixed_temps = [2.0, 2.2, 3.0, 4.0] #default fixed temperatures to analyze

#Information on which functions to fit to the different RF curves
fit_funs_bt = ['BCS','BCS','BCS','p2','p2','p2','p2','p2','p2','p2']
fit_funs_at = ['BCS','BCS','BCS','p4','p3','p2','BCS','p3','p4','p4']

skiplines = 0 #The default number of lines to skip at the beginning of the table is zero

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

#Subtract 1 because Python indexing starts at 0
Qcol -= 1
Tcol -= 1
Eacccol -= 1
cavity -= 1

#Get frequency and G factor values
freq = frequencies[cavity]
G = G_vals[cavity]

#Create empty lists to hold the data
Qdata = []
Tdata = []
Eaccdata = []
Rssdata = []
"""---------------------------------------------------------------------------------------
Read in the data from the input file
---------------------------------------------------------------------------------------"""
#Read the input file line by line and extract data
line_cnt = 1 #line count, it will start with line 1
for line in file:

    if line_cnt <= skiplines: #skip the indecated number of lines at the beginning
        line_cnt += 1
        continue

    line = line.strip()
    columns = line.split()

    if columns == []: #skip any columns that are completely empty
        line_cnt += 1
        continue

    try: #add data to the lists
        if float(columns[Qcol]) > 0:
            Qdata.append(float(columns[Qcol]))
        else:
            print("Warning: Qo value is less than 0 on line", line_cnt, "skipping line", line_cnt)
            line_cnt += 1
            continue
    except IndexError:
        print("Warning: index error on line", line_cnt, "in the Q column. Skipping line", line_cnt)
        line_cnt += 1
        continue
    except ValueError:
        print("Value on line", line_cnt, "in the Q column is not a number. Skipping line", line_cnt)
        line_cnt += 1
        continue

    try: #Skip an entire line if any of the columns has an error so that the lists are of the same length
        Tdata.append(float(columns[Tcol]))
    except IndexError:
        print("Warning: index error on line", line_cnt, "in the Temperature column. Skipping line", line_cnt)
        Qdata.pop()
        line_cnt += 1
        continue
    except ValueError:
        print("Value on line", line_cnt, "in the Temperature column is not a number. Skipping line", line_cnt)
        Qdata.pop()
        line_cnt += 1
        continue

    try:
        Eaccdata.append(float(columns[Eacccol]))
    except IndexError:
        print("Warning: index error on line", line_cnt, "in the Eacc column. Skipping line", line_cnt)
        Qdata.pop()
        Tdata.pop()
        line_cnt += 1
        continue
    except ValueError:
        print("Value on line", line_cnt, "in the Eacc column is not a number. Skipping line", line_cnt)
        Qdata.pop()
        Tdata.pop()
        line_cnt += 1
        continue

    line_cnt += 1

file.close()
#print(len(Qdata), len(Tdata), len(Eaccdata))

"""---------------------------------------------------------------------------------------
For each accelerating field ramp up, calculate the feild corrected Rs values
---------------------------------------------------------------------------------------"""
#Get the non-corrected Rs* data
for Q in Qdata:
    Rssdata.append(G/Q)

ramp_Eacc = []
ramp_Rss = []
Rs = [] #Array for corrected Rs data

#Function for doing the Rs* to Rs correction using the beta factors
def Rs_correction(Bp, Rparam):
    return beta3_all[cavity]*Rparam[0]*Bp**3 + beta2_all[cavity]*Rparam[1]*Bp**2 + beta1_all[cavity]*Rparam[2]*Bp + beta0_all[cavity]*Rparam[3]

#This loop gets the corrected Rs data and puts it in Rs
for i in range(0,len(Eaccdata)):

    if i == 0:
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])
        continue

    if i == (len(Eaccdata)-1): #if you have reached the last data point

        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])

        #Do the fit here
        coef = np.polyfit(ramp_Eacc,ramp_Rss,3) #3rd degree polynomial
        #print(coef)
        #fit = np.poly1d(coef)
        #print(fit)

        for j in range(0,len(ramp_Rss)):
            Rs.append(Rs_correction(ramp_Eacc[j], coef))
        #print("last ramp", ramp_Eacc[0])

    if Eaccdata[i-1] < Eaccdata[i]:
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])
    else:

        #Do the fit here
        coef = np.polyfit(ramp_Eacc,ramp_Rss,3) #3rd degree polynomial
        #print(coef)
        #fit = np.poly1d(coef)
        #print(fit)

        for j in range(0,len(ramp_Rss)):
            Rs.append(Rs_correction(ramp_Eacc[j], coef))
        #print("some ramp", ramp_Eacc[0])

        #Get ready to start the next ramp up fit
        ramp_Eacc.clear()
        ramp_Rss.clear()
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])

#print(ramp_Eacc)
#print(len(Tdata), len(Rssdata), len(Eaccdata), len(Rs))
"""---------------------------------------------------------------------------------------
Separate the data before and after the superfluid helium transition
---------------------------------------------------------------------------------------"""
#The _bt means before transition, as in before the superfluid transition of helium
Qdata_bt = []
Tdata_bt = []
Eaccdata_bt = []
Rs_bt = [] #Rs is corrected surface resistance data

#The _at means after transition, as in after the superfluid transition of helium
Qdata_at = []
Tdata_at = []
Eaccdata_at = []
Rs_at = []

for i in range(0, len(Tdata)):

    if Tdata[i] >= 2.174:
        Qdata_bt.append(Qdata[i])
        Tdata_bt.append(Tdata[i])
        Eaccdata_bt.append(Eaccdata[i])
        Rs_bt.append(Rs[i])

    elif Tdata[i] < 2.174:
        Qdata_at.append(Qdata[i])
        Tdata_at.append(Tdata[i])
        Eaccdata_at.append(Eaccdata[i])
        Rs_at.append(Rs[i])
"""---------------------------------------------------------------------------------------
Separate the data by different field amplitude values
---------------------------------------------------------------------------------------"""
Inv_Tdata_sep_bt = [] #list of lists for inverse temperature data separated by field amplitude before superfluid helium transition (bt)
Tdata_sep_bt = [] #list of lists for temperature data separated by field amplitude before superfluid helium transition (bt)
Rs_sep_bt = [] #list of lists for surface resistance data in nano-ohms separated by field amplitude before superfluid helium transition (bt)

Inv_Tdata_sep_at = [] #list of lists for inverse temperature data separated by field amplitude after superfluid helium transition (at)
Tdata_sep_at = [] #list of lists for temperature data separated by field amplitude after superfluid helium transition (at)
Rs_sep_at = [] #list of lists for surface resistance data in nano-ohms separated by field amplitude after superfluid helium transition (at)

for value in FieldValues:

    Inv_Tdata_sep_bt.append([])
    Tdata_sep_bt.append([])
    Rs_sep_bt.append([])

    Inv_Tdata_sep_at.append([])
    Tdata_sep_at.append([])
    Rs_sep_at.append([])

#Make two lists containing lists of the inverse temperature data and the
#corrected surface resistance data. Each sub-list corresponds to a different
#field amplitude. Do this once for data before superfluid helium transition,
#and again for after superfluid helium transition
for i in range(0,len(Eaccdata_bt)):

    for j in range(0,len(FieldValues)):

        if (Eaccdata_bt[i] < FieldValues[j]+0.5) and (Eaccdata_bt[i] >= FieldValues[j]-0.5):
            Inv_Tdata_sep_bt[j].append(1/Tdata_bt[i])
            Tdata_sep_bt[j].append(Tdata_bt[i])
            Rs_sep_bt[j].append(Rs_bt[i]*(10**9))

for i in range(0,len(Eaccdata_at)):

    for j in range(0,len(FieldValues)):

        if (Eaccdata_at[i] < FieldValues[j]+0.5) and (Eaccdata_at[i] >= FieldValues[j]-0.5):
            Inv_Tdata_sep_at[j].append(1/Tdata_at[i])
            Tdata_sep_at[j].append(Tdata_at[i])
            Rs_sep_at[j].append(Rs_at[i]*(10**9))

"""---------------------------------------------------------------------------------------
Get values for the superfluid transition of helium
---------------------------------------------------------------------------------------"""
"""
SF_field_vals = [] #list of rf field values for ramp before and after sf helium transition
SF_Rs_vals = [] #Corresponding rf field values
SF_T_vals = []

ramp_vals_cnt_bt = 0
ramp_vals_cnt = 0

for i in range(0, len(Tdata)):

    if i==0:
        ramp_vals_cnt += 1
        continue

    if Eaccdata[i-1] < Eaccdata[i]:
        ramp_vals_cnt += 1
        continue

    elif (Eaccdata[i-1] > Eaccdata[i]) and Tdata[i] > 2.175:
        ramp_vals_cnt_bt = ramp_vals_cnt
        ramp_vals_cnt = 1
        continue

    elif (Eaccdata[i-1] > Eaccdata[i]) and Tdata[i] <= 2.175:

        for j in range((i-ramp_vals_cnt-ramp_vals_cnt_bt), i):
            SF_field_vals.append(Eaccdata[j])
            SF_Rs_vals.append(Rs[j])
            SF_T_vals.append(Tdata[j])

        n = i
        while Eaccdata[n] < Eaccdata[n+1]:
            SF_field_vals.append(Eaccdata[n])
            SF_Rs_vals.append(Rs[n])
            SF_T_vals.append(Tdata[n])
            n += 1
        SF_field_vals.append(Eaccdata[n])
        SF_Rs_vals.append(Rs[n])
        SF_T_vals.append(Tdata[n])
        break

file_sf = open('sf_data.csv', 'a')
file_sf.write(file.name)
file_sf.write('\n')

for k in range(0, len(SF_field_vals)):
    file_sf.write(str(SF_field_vals[k]))
    file_sf.write(',')
    file_sf.write(str(SF_Rs_vals[k]))
    file_sf.write(',')
    file_sf.write(str(SF_T_vals[k]))
    file_sf.write('\n')

file_sf.close()
"""
"""---------------------------------------------------------------------------------------
Fit the data to the RBCS formula or a polynomial and plot the results
---------------------------------------------------------------------------------------"""

#RBCS fit function
def BCS(T, a0, a1, Rres, f=freq, Tc=9.25):
    kB = 8.617333262145*10**(-5) #Boltzman constant
    hbar = 6.582119569*10**(-16)

    a1T = a1*np.sqrt(np.cos((np.pi/2)*(T/Tc)**2))
    C = 8/(np.exp(0.5772156649))

    return (10**9)*(a0/T)*np.log(C*kB*T/(2*np.pi*hbar*f*10**6))*np.exp(-a1T*Tc/T) + Rres

#Second order polynomial fit function
def Poly2(T, a, b, c):
    return (a*T*T + b*T + c)

#Third order polynomial fit function
def Poly3(T, a, b, c, d):
    return (a*T*T*T + b*T*T + c*T + d)

#Fourth order polynomial fit function
def Poly4(T, a, b, c, d, e):
    return (a*T*T*T*T + b*T*T*T + c*T*T + d*T + e)

#Make the models
fmodel = Model(BCS)
p2model = Model(Poly2)
p3model = Model(Poly3)
p4model = Model(Poly4)

fig1, ax1 = plt.subplots(nrows=1, ncols=1)
fig2, ax2 = plt.subplots(nrows=1, ncols=1)
fig3, ax3 = plt.subplots(nrows=1, ncols=1)

colors = ['b','orange', 'g', 'r', 'c', 'm', 'y', 'salmon', 'brown', 'lawngreen' , '0.4', '0.8' ]
shapes = ['^', 's', 'P', '*', '+', 'd', 'x']

legend_entries = []
temps_legend = []

#Values of the total surface resistance for fixed temperatures calculated from the fits
Rs_fixed_temps = []

#Make a list of lists to store the total surface resistance for fixed temperature values
#each sub-list corresponds to a different fixed temperature
for i in range(0,len(fixed_temps)):
    Rs_fixed_temps.append([])

for i in range(0,len(FieldValues)):

    name = str(FieldValues[i]) + " mT"
    legend_entries.append(name)

for i in range(0, len(fixed_temps)):
    temp_str = str(fixed_temps[i])
    temps_legend.append(temp_str)

if len(legend_entries) > (len(colors)+len(shapes)):
    print("Warning: This program wasn't expecting more than 19 different field amplitudes. Some field amplitude data sets will be indistiguishable on the plot, as they will be plotted as black dots.")

for i in range(0,len(legend_entries)):

    #Make inital guesses for the residual resistance
    res_min_bt = min(Rs_sep_bt[i])
    res_min_at = min(Rs_sep_at[i])

    if fit_funs_bt[i] == 'BCS':
        #Fit each feild amplitude data set to the RBCS formula
        params_bt = fmodel.make_params(a0=0.001, a1=1.5, Rres=res_min_bt, f=freq, Tc=9.25)
        params_bt['f'].vary = False
        params_bt['Tc'].vary = False

        #fit is made in this line
        result_bt = fmodel.fit(Rs_sep_bt[i], params_bt, T=Tdata_sep_bt[i])

        #Get the parameters calculated from the fits
        a0_fit_bt = result_bt.best_values['a0']
        a1_fit_bt = result_bt.best_values['a1']
        Rres_fit_bt = result_bt.best_values['Rres']

    elif fit_funs_bt[i] == 'p2':
        params_bt = p2model.make_params(a=2, b=1, c=res_min_bt)
        result_bt = p2model.fit(Rs_sep_bt[i], params_bt, T=Tdata_sep_bt[i])

        #Get the parameters calculated from the fits
        a_fit_bt = result_bt.best_values['a']
        b_fit_bt = result_bt.best_values['b']
        c_fit_bt = result_bt.best_values['c']

    elif fit_funs_bt[i] == 'p3':
        params_bt = p3model.make_params(a=3, b=2, c=1, d=res_min_bt)
        result_bt = p3model.fit(Rs_sep_bt[i], params_bt, T=Tdata_sep_bt[i])

        #Get the parameters calculated from the fits
        a_fit_bt = result_bt.best_values['a']
        b_fit_bt = result_bt.best_values['b']
        c_fit_bt = result_bt.best_values['c']
        d_fit_bt = result_bt.best_values['d']

    elif fit_funs_bt[i] == 'p4':
        params_bt = p4model.make_params(a=4, b=3, c=2, d=1, e=res_min_bt)
        result_bt = p4model.fit(Rs_sep_bt[i], params_bt, T=Tdata_sep_bt[i])

        #Get the parameters calculated from the fits
        a_fit_bt = result_bt.best_values['a']
        b_fit_bt = result_bt.best_values['b']
        c_fit_bt = result_bt.best_values['c']
        d_fit_bt = result_bt.best_values['d']
        e_fit_bt = result_bt.best_values['e']

    if fit_funs_at[i] == 'BCS':
        #Fit each feild amplitude data set to the RBCS formula
        params_at = fmodel.make_params(a0=0.001, a1=1.5, Rres=res_min_at, f=freq, Tc=9.25)
        params_at['f'].vary = False
        params_at['Tc'].vary = False

        #fit is made in this line
        result_at = fmodel.fit(Rs_sep_at[i], params_at, T=Tdata_sep_at[i])

        #Get the parameters calculated from the fits
        a0_fit_at = result_at.best_values['a0']
        a1_fit_at = result_at.best_values['a1']
        Rres_fit_at = result_at.best_values['Rres']

    elif fit_funs_at[i] == 'p2':
        params_at = p2model.make_params(a=2, b=1, c=res_min_at)
        result_at = p2model.fit(Rs_sep_at[i], params_at, T=Tdata_sep_at[i])

        #Get the parameters calculated from the fits
        a_fit_at = result_at.best_values['a']
        b_fit_at = result_at.best_values['b']
        c_fit_at = result_at.best_values['c']

    elif fit_funs_at[i] == 'p3':
        params_at = p3model.make_params(a=3, b=2, c=1, d=res_min_at)
        result_at = p3model.fit(Rs_sep_at[i], params_at, T=Tdata_sep_at[i])

        #Get the parameters calculated from the fits
        a_fit_at = result_at.best_values['a']
        b_fit_at = result_at.best_values['b']
        c_fit_at = result_at.best_values['c']
        d_fit_at = result_at.best_values['d']

    elif fit_funs_at[i] == 'p4':
        params_at = p4model.make_params(a=4, b=3, c=2, d=1, e=res_min_at)
        result_at = p4model.fit(Rs_sep_at[i], params_at, T=Tdata_sep_at[i])

        #Get the parameters calculated from the fits
        a_fit_at = result_at.best_values['a']
        b_fit_at = result_at.best_values['b']
        c_fit_at = result_at.best_values['c']
        d_fit_at = result_at.best_values['d']
        e_fit_at = result_at.best_values['e']

    #Get the data for the total surface resistance calculated from the fits for
    #before and after the superfluid transition
    for j in range(0, len(fixed_temps)):

        if fixed_temps[j] >= 2.175:
            if fit_funs_bt[i] == 'BCS':
                Rs_fixed_temps[j].append(BCS(fixed_temps[j],a0_fit_bt,a1_fit_bt,Rres_fit_bt))

            elif fit_funs_bt[i] == 'p2':
                Rs_fixed_temps[j].append(Poly2(fixed_temps[j], a_fit_bt, b_fit_bt, c_fit_bt))

            elif fit_funs_bt[i] == 'p3':
                Rs_fixed_temps[j].append(Poly3(fixed_temps[j], a_fit_bt, b_fit_bt, c_fit_bt, d_fit_bt))

            elif fit_funs_bt[i] == 'p4':
                Rs_fixed_temps[j].append(Poly4(fixed_temps[j], a_fit_bt, b_fit_bt, c_fit_bt, d_fit_bt, e_fit_bt))

        elif fixed_temps[j] < 2.175:
            if fit_funs_at[i] == 'BCS':
                Rs_fixed_temps[j].append(BCS(fixed_temps[j],a0_fit_at,a1_fit_at,Rres_fit_at))

            elif fit_funs_at[i] == 'p2':
                Rs_fixed_temps[j].append(Poly2(fixed_temps[j], a_fit_at, b_fit_at, c_fit_at))

            elif fit_funs_at[i] == 'p3':
                Rs_fixed_temps[j].append(Poly3(fixed_temps[j], a_fit_at, b_fit_at, c_fit_at, d_fit_at))

            elif fit_funs_at[i] == 'p4':
                Rs_fixed_temps[j].append(Poly4(fixed_temps[j], a_fit_at, b_fit_at, c_fit_at, d_fit_at, e_fit_at))

    #plot the data points and fit lines
    if i < len(colors):
        ax1.plot(Inv_Tdata_sep_bt[i],Rs_sep_bt[i], marker='o', linestyle='none', markersize=4, color=colors[i], label=legend_entries[i])
        ax1.plot(Inv_Tdata_sep_at[i],Rs_sep_at[i], marker='o', linestyle='none', markersize=4, color=colors[i])

        ax1.plot(Inv_Tdata_sep_bt[i], result_bt.best_fit, marker='None', linestyle='--',color=colors[i])
        ax1.plot(Inv_Tdata_sep_at[i], result_at.best_fit, marker='None', linestyle='--',color=colors[i])

        ax2.plot(Inv_Tdata_sep_bt[i],((result_bt.residual))*100/Rs_sep_bt[i], marker='o', markersize=3, color=colors[i], label=legend_entries[i])
        ax2.plot(Inv_Tdata_sep_at[i],((result_at.residual))*100/Rs_sep_at[i], marker='o', markersize=3, color=colors[i])

    elif i < (len(colors)+len(shapes)):
        ax1.plot(Inv_Tdata_sep[i],Rs_sep[i], marker=shapes[i-len(colors)], linestyle='none', markersize=4, color='black', label=legend_entries[i])
        ax1.plot(Inv_Tdata_sep[i], result.best_fit, marker='None', linestyle='--',color='black')
        ax2.plot(Inv_Tdata_sep[i],((result.residual))*100/Rs_sep[i], marker=shapes[i-len(colors)], markersize=3, color='black', label=legend_entries[i])
    else:
        ax1.plot(Inv_Tdata_sep[i],Rs_sep[i], marker='o', linestyle='none', markersize=4, color='black', label=legend_entries[i])
        ax1.plot(Inv_Tdata_sep[i], result.best_fit, marker='None', linestyle='--',color='black')
        ax2.plot(Inv_Tdata_sep[i],((result.residual))*100/Rs_sep[i], marker='o', markersize=3, color='black', label=legend_entries[i])

for i in range(0, len(fixed_temps)):
    ax3.plot(FieldValues, Rs_fixed_temps[i], marker='o', linestyle='none', markersize=4, label=temps_legend[i])

x_formatter = FixedFormatter([r'10$^{-1}$', r'9$^{-1}$', r'8$^{-1}$', r'7$^{-1}$', r'6$^{-1}$', r'5$^{-1}$',
 r'4.5$^{-1}$', r'4.0$^{-1}$', r'3.5$^{-1}$', r'3.0$^{-1}$', r'2.5$^{-1}$', r'2.2$^{-1}$', r'2.0$^{-1}$', r'1.8$^{-1}$', r'1.7$^{-1}$'])
x_locator = FixedLocator([0.1, 1/9, 0.125, 1/7, 1/6, 0.2, 1/4.5, 0.25, 1/3.5, 1/3, 0.4, 1/2.2, 0.5, 1/1.8, 1/1.7])

ax1.set_yscale('log')
ax1.set_xlabel(r'Inverse Temparature [K$^{-1}$]')
ax1.set_ylabel(r'R$_s[n\Omega]$')
ax1.xaxis.set_major_formatter(x_formatter)
ax1.xaxis.set_major_locator(x_locator)
ax1.legend(title='Field Amplitude')
ax1.set_title("Surface Resistance vs Inverse Temperature")
ax1.grid(True)

ax2.set_xlabel(r'Inverse Temparature [K$^{-1}$]')
ax2.set_ylabel(r"Residual/R$_s$ [%]")
ax2.xaxis.set_major_formatter(x_formatter)
ax2.xaxis.set_major_locator(x_locator)
ax2.set_title(r"Percent Difference between observed R$_s$ and values calculated from fit parameters")
ax2.legend(title='Field Amplitude')
ax2.grid(True)

#ax3.set_yscale('log')
ax3.set_xlabel('RF field [mT]')
ax3.set_ylabel(r'R$_s[n\Omega]$')
ax3.legend(title='Temperature [K]')
ax3.set_title('Total Surface Resistance vs RF Field for Fixed Temperatures')
ax3.grid(True)

plt.show()
