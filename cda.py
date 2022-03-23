"""
 *========================================================================================
 *
 *       File name:  cda.py
 *
 *    Description:  Takes Q vs T data collected during a cooldown and extracts
 *                  Rbcs and R0 as a function of field amplitude
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
import pandas as pd
from lmfit import Model

"""---------------------------------------------------------------------------------------
Get all the input parameters ready
---------------------------------------------------------------------------------------"""
#These are the beta values to use in the fit function corrections when converting Rs* to Rs
#The six entries in the lists coR0pond to the six different coaxial cavities
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

fixed_temps = [1.9, 2.2, 2.5, 2.8, 3.1 ,3.4, 3.7, 4.0] #default fixed temperatures to analyze

#Information on which functions to fit to the different RF curves
fit_funs = ['BCS','BCS','BCS','BCS','BCS','BCS','BCS','BCS','BCS','BCS','BCS']

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
SWR = 1.080124

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
    Rssdata.append((G/Qdata[i])*10**9)
    Bp_err.append(Eaccdata[i]*(SWR-1)/4)
    Rss_err.append(Rssdata[i]*(SWR-1)/2)
    weights.append(1 - (Rss_err[i]/Rssdata[i]))

ramp_Eacc = []
ramp_Rss = []
ramp_weights = []

Rs = [] #Array for corrected Rs data
Rs_err = [] #Array for uncertainties in corrected Rs data

#Function for doing the Rs* to Rs correction using the beta factors
def Rs_correction(Bp, Rparam):
    return beta3_all[cavity]*Rparam[0]*Bp**3 + beta2_all[cavity]*Rparam[1]*Bp**2 + beta1_all[cavity]*Rparam[2]*Bp + beta0_all[cavity]*Rparam[3]

#Do a reverse Rs to Rs* correction
def Rev_Rs_correction(Rs_ramp,Bp):

    Bp_ramp = FieldValues
    #print(Bp_ramp)
    #print(Rs_ramp)
    coef = np.polyfit(Bp_ramp,Rs_ramp,3) #3rd degree polynomial

    Rparam = coef

    Rss = (Rparam[0]*Bp**3)/(beta3_all[cavity]) + (Rparam[1]*Bp**2)/(beta2_all[cavity]) + (Rparam[2]*Bp)/(beta1_all[cavity]) + (Rparam[3])/(beta0_all[cavity])

    return Rss

def Delta_Rs(Bp,Rparam,delta_Rparam,delta_Bp):

    dRs_squared = 0
    beta_vals = [beta3_all[cavity], beta2_all[cavity], beta1_all[cavity], beta0_all[cavity]]
    alpha_i_powers = [3,2,1,0]

    for i in range(0,len(Rparam)):

        Rs_val = beta_vals[i]*Rparam[i]*Bp**(alpha_i_powers[i])
        sqrt_val = (delta_Rparam[i]/Rparam[i])**2 + (alpha_i_powers[i]*delta_Bp/Bp)**2

        dRs_squared += (Rs_val*np.sqrt(sqrt_val))**2

    dRs = np.sqrt(dRs_squared)

    return dRs

#This loop gets the corrected Rs data and puts it in Rs
for i in range(0,len(Eaccdata)):

    if i == 0:
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])
        ramp_weights.append(weights[i])
        continue

    if i == (len(Eaccdata)-1): #if you have reached the last data point

        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])
        ramp_weights.append(weights[i])

        #Do the fit here
        coef,cov = np.polyfit(ramp_Eacc,ramp_Rss,3, w=ramp_weights, cov=True) #3rd degree polynomial
        #print(coef)
        #fit = np.poly1d(coef)
        #print(fit)

        for j in range(0,len(ramp_Rss)):
            Rs.append(Rs_correction(ramp_Eacc[j], coef))

            Rparam_err = np.sqrt(np.diag(cov))

            Bp_err = ramp_Eacc[j]*(SWR-1)/4

            Rs_err.append(Delta_Rs(ramp_Eacc[j], coef, Rparam_err, Bp_err))
        #print("last ramp", ramp_Eacc[0])

    if Eaccdata[i-1] < (Eaccdata[i] + 0.5):
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])
        ramp_weights.append(weights[i])
    else:

        #Do the fit here
        #print(len(ramp_Eacc),Eaccdata[i])
        try:
            #print(len(ramp_Eacc), len(ramp_Rss))
            coef,cov = np.polyfit(ramp_Eacc,ramp_Rss,3, w=ramp_weights, cov=True) #3rd degree polynomial
            #print(coef)

            #fit = np.poly1d(coef)
            #print(fit)
        except ValueError:
            ramp_Eacc.append(Eaccdata[i])
            ramp_Rss.append(Rssdata[i])
            ramp_weights.append(weights[i])
            continue

        for j in range(0,len(ramp_Rss)):
            Rs.append(Rs_correction(ramp_Eacc[j], coef))

            Rparam_err = np.sqrt(np.diag(cov))

            Bp_err = ramp_Eacc[j]*(SWR-1)/4

            Rs_err.append(Delta_Rs(ramp_Eacc[j], coef, Rparam_err, Bp_err))
        #print("some ramp", ramp_Eacc[0])

        #Get ready to start the next ramp up fit
        ramp_Eacc.clear()
        ramp_Rss.clear()
        ramp_weights.clear()
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])
        ramp_weights.append(weights[i])

#print(Rs_err)
#print(ramp_Eacc)
#print(len(Tdata), len(Rssdata), len(Eaccdata), len(Rs), len(Rs_err))
"""---------------------------------------------------------------------------------------
Separate the data by different field amplitude values
---------------------------------------------------------------------------------------"""

Inv_Tdata_sep = [] #list of lists for inverse temperature data separated by field amplitude
Tdata_sep = [] #list of listst for temperature data separated by field amplitude
Rs_sep = [] #list of lists for surface resistance data in nano-ohms separated by field amplitude
Rs_err_sep = []
Rs_sep_log = [] #list of lists for the logarithm of the surface resistance data in nano-ohms separated by field amplitude
Rs_sep_ln = [] #list of lists for the natural logarithm of the surface resistance data in nano-ohms separated by field amplitude
weights_sep = []

for value in FieldValues:
    Inv_Tdata_sep.append([])
    Tdata_sep.append([])
    Rs_sep.append([])
    Rs_err_sep.append([])
    Rs_sep_log.append([])
    Rs_sep_ln.append([])
    weights_sep.append([])

#Make two lists containing lists of the inverse temperature data and the
#corrected surface resistance data. Each sub-list coR0ponds to a different
#field amplitude

for i in range(0,len(Eaccdata)):

    for j in range(0,len(FieldValues)):

        if (Eaccdata[i] < FieldValues[j]+0.5) and (Eaccdata[i] > FieldValues[j]-0.5):
            Inv_Tdata_sep[j].append(1/Tdata[i])
            Tdata_sep[j].append(Tdata[i])
            Rs_sep[j].append(Rs[i])
            Rs_err_sep[j].append(Rs_err[i])
            Rs_sep_log[j].append(np.log10(Rs[i]))
            Rs_sep_ln[j].append(np.log(Rs[i]))
            weights_sep[j].append(weights[i])

"""---------------------------------------------------------------------------------------
Fit the data to the RBCS formula and plot the results
---------------------------------------------------------------------------------------"""
#RBCS fit function
def BCS(T, a0, a1, R0, DeltaRs, f=freq, Tc=9.25):
    kB = 8.617333262145*10**(-5) #Boltzman constant
    hbar = 6.582119569*10**(-16)
    Step = np.sign(T-2.175)

    try:
        for i in range(0,len(Step)):
            if Step[i] == -1:
                Step[i] = 1
            elif Step[i] == 1:
                Step[i] = 0

    except TypeError:
        if Step == -1:
            Step = 1
        elif Step == 1:
            Step = 0

    a1T = a1*np.sqrt(np.cos((np.pi/2)*(T/Tc)**2))
    C = 8/(np.exp(0.5772156649))

    result = ((10**9)*(a0/T)*np.log(C*kB*T/(2*np.pi*hbar*f*10**6))*np.exp(-a1T*Tc/T) + R0 + DeltaRs*Step)
    #print('func', result, a0, a1, R0, DeltaRs)

    return result
    #return np.log(result)

def Plotting_BCS_bt(T, a0, a1, R0, DeltaRs, f=freq, Tc=9.25):
    kB = 8.617333262145*10**(-5) #Boltzman constant
    hbar = 6.582119569*10**(-16)

    a1T = a1*np.sqrt(np.cos((np.pi/2)*(T/Tc)**2))
    C = 8/(np.exp(0.5772156649))

    result = ((10**9)*(a0/T)*np.log(C*kB*T/(2*np.pi*hbar*f*10**6))*np.exp(-a1T*Tc/T) + R0)
    #print(result, a0, a1, R0, DeltaRs)

    return result

def Plotting_BCS_at(T, a0, a1, R0, DeltaRs, f=freq, Tc=9.25):
    kB = 8.617333262145*10**(-5) #Boltzman constant
    hbar = 6.582119569*10**(-16)

    a1T = a1*np.sqrt(np.cos((np.pi/2)*(T/Tc)**2))
    C = 8/(np.exp(0.5772156649))

    result = ((10**9)*(a0/T)*np.log(C*kB*T/(2*np.pi*hbar*f*10**6))*np.exp(-a1T*Tc/T) + R0 + DeltaRs)
    #print(result, a0, a1, R0, DeltaRs)

    return result

#Second order polynomial fit function
def Poly2(T, a, b, c):
    return (a*T*T + b*T + c)

#Third order polynomial fit function
def Poly3(T, a, b, c, d):
    return (a*T*T*T + b*T*T + c*T + d)

#Fourth order polynomial fit function
def Poly4(T, a, b, c, d, e):
    return (a*T*T*T*T + b*T*T*T + c*T*T + d*T + e)

def fixed_T_Poly(B, alpha, beta, gamma):
    return (alpha + beta*B + gamma*B*B)

def fixed_T_quad(B, R0q, gammaq):
    return (R0q + R0q*gammaq*(B/100)*(B/100))

#Make the models
fmodel = Model(BCS)
p2model = Model(Poly2)
p3model = Model(Poly3)
p4model = Model(Poly4)
fixed_T_model = Model(fixed_T_Poly)
fixed_T_quad_model = Model(fixed_T_quad)

fig1, ax1 = plt.subplots(nrows=1, ncols=1)
fig2, ax2 = plt.subplots(nrows=1, ncols=1)
fig3, ax3 = plt.subplots(nrows=1, ncols=1)
fig4, ax4 = plt.subplots(nrows=1, ncols=1)

colors = ['b','orange', 'g', 'r', 'c', 'm', 'y', 'salmon', 'brown', 'lawngreen' , '0.4', '0.8' ]
shapes = ['^', 's', 'P', '*', '+', 'd', 'x']

legend_entries = []
temps_legend = []

#Values of the total surface resistance for fixed temperatures calculated from the fits
Rs_fixed_temps = []

#Values of the non-corrected surface resistance for fixed temperatures calculated from the fits
Rss_fixed_temps = []
Q_fixed_temps = []

#Temperature values for the fit plotting
T_vals = np.linspace(1.9, 4.5, 50)
Inv_T_vals = np.linspace(1/1.9, 1/4.5, 50)
Rs_for_T_vals = []
Rss_for_T_vals = []

#Make a list of lists to store the total surface resistance for fixed temperature values
#each sub-list coR0ponds to a different fixed temperature
for i in range(0,len(fixed_temps)):
    Rs_fixed_temps.append([])
    Rss_fixed_temps.append([])
    Q_fixed_temps.append([])

for i in range(0,len(FieldValues)):

    name = str(FieldValues[i]) + " mT"
    legend_entries.append(name)

    Rs_for_T_vals.append([])
    Rss_for_T_vals.append([])

for i in range(0, len(fixed_temps)):
    temp_str = str(fixed_temps[i])
    temps_legend.append(temp_str)

if len(legend_entries) > (len(colors)+len(shapes)):
    print("Warning: This program wasn't expecting more than 19 different field amplitudes. Some field amplitude data sets will be indistiguishable on the plot, as they will be plotted as black dots.")


file_fit_params = open('fit_params.txt', 'w')
#file_fit_params.write("a0 , a1 , R0")
#file_fit_params.write('\n')

for i in range(0,len(legend_entries)):

    #Make inital guesses for the residual resistance and Delta Rs
    res_min_guess = min(Rs_sep[i])
    #print('Res min: ', res_min_guess)
    #print()
    DeltaRs_guess = i/10

    if fit_funs[i] == 'BCS':
        #Fit each feild amplitude data set to the RBCS formula
        params = fmodel.make_params(a0=1, a1=1.5, R0=res_min_guess, DeltaRs=DeltaRs_guess, f=freq, Tc=9.25)
        params['f'].vary = False
        params['Tc'].vary = False
        #params['DeltaRs'].set(min=-20, max=0)
        params['a0'].set(min=0)
        #params['a1'].set(min=0,max=10)
        params['R0'].set(min=(res_min_guess-2),max=(res_min_guess))
        #if i == 1:
            #params['DeltaRs'].set(min=-2, max=1)
            #params['a1'].set(max=2)
        #if i==3:
            #params['DeltaRs'].set(min=-4,max=1)
        #if i==5:
            #params['DeltaRs'].set(min=-7,max=1)

        #fit is made in this line
        result = fmodel.fit(Rs_sep[i], params, T=Tdata_sep[i])
        #print(Tdata_sep[i])
        #print(Rs_sep_ln[i])
        #print(result.best_fit)

        #Get the parameters calculated from the fits
        a0_fit = result.best_values['a0']
        a1_fit = result.best_values['a1']
        R0_fit = result.best_values['R0']
        DeltaRs_fit = result.best_values['DeltaRs']
        #print(a0_fit, a1_fit, R0_fit, DeltaRs_fit)
        #print('delta Rs: ', round(np.exp(-1*DeltaRs_fit),4))
        #print(round((-1*DeltaRs_fit),2), '+/-', round(result.params['DeltaRs'].stderr,2))
        #print(round((-1*DeltaRs_fit),2))
        #print(round(result.params['DeltaRs'].stderr,2))
        #print(R0_fit)

        #print(result.fit_report())

        #file_fit_params.write(str(a0_fit))
        #file_fit_params.write(',')
        #file_fit_params.write(str(a1_fit))
        #file_fit_params.write(',')
        #file_fit_params.write(str(R0_fit))
        #file_fit_params.write('\n')

        file_fit_params.write('Bp: ')
        file_fit_params.write(str(legend_entries[i]))
        file_fit_params.write(" [mT]")
        file_fit_params.write('\n')
        file_fit_params.write(result.fit_report())
        file_fit_params.write('\n')
        file_fit_params.write('\n')

    elif fit_funs[i] == 'p2':
        params = p2model.make_params(a=2, b=1, c=res_min)
        result = p2model.fit(Rs_sep[i], params, T=Tdata_sep[i])

        #Get the parameters calculated from the fits
        a_fit = result.best_values['a']
        b_fit = result.best_values['b']
        c_fit = result.best_values['c']

    elif fit_funs[i] == 'p3':
        params = p3model.make_params(a=3, b=2, c=1, d=res_min)
        result = p3model.fit(Rs_sep[i], params, T=Tdata_sep[i])

        #Get the parameters calculated from the fits
        a_fit = result.best_values['a']
        b_fit = result.best_values['b']
        c_fit = result.best_values['c']
        d_fit = result.best_values['d']

    elif fit_funs[i] == 'p4':
        params = p4model.make_params(a=4, b=3, c=2, d=1, e=res_min)
        result = p4model.fit(Rs_sep[i], params, T=Tdata_sep[i])

        #Get the parameters calculated from the fits
        a_fit = result.best_values['a']
        b_fit = result.best_values['b']
        c_fit = result.best_values['c']
        d_fit = result.best_values['d']
        e_fit = result.best_values['e']

    #Get the data for the total surface resistance calculated from the fits
    for j in range(0, len(fixed_temps)):

        if fit_funs[i] == 'BCS':
            Rs_fixed_temps[j].append(BCS(fixed_temps[j],a0_fit,a1_fit,R0_fit,DeltaRs_fit))
            #print(round(R0_fit,4))

        #elif fit_funs[i] == 'p2':
            #Rs_fixed_temps[j].append(Poly2(fixed_temps[j], a_fit, b_fit, c_fit))

        #elif fit_funs[i] == 'p3':
            #Rs_fixed_temps[j].append(Poly3(fixed_temps[j], a_fit, b_fit, c_fit, d_fit))

        #elif fit_funs[i] == 'p4':
            #Rs_fixed_temps[j].append(Poly4(fixed_temps[j], a_fit, b_fit, c_fit, d_fit, e_fit))

    #for k in range(0, len(T_vals)):

        #Rs_for_T_vals[i].append(BCS(T_vals[k], a0_fit,a1_fit,R0_fit,DeltaRs_fit))

    #plot the data points and fit lines
    if i < len(colors):
        Inv_Tdata_sep_bt = []
        Tdata_sep_bt = []
        Inv_Tdata_sep_at = []
        Tdata_sep_at = []
        for j in range(0,len(Inv_Tdata_sep[i])):
            if Inv_Tdata_sep[i][j] < (1/2.175):
                Inv_Tdata_sep_bt.append(Inv_Tdata_sep[i][j])
                Tdata_sep_bt.append(Tdata_sep[i][j])
            else:
                Inv_Tdata_sep_at.append(Inv_Tdata_sep[i][j])
                Tdata_sep_at.append(Tdata_sep[i][j])

        x_bt = np.linspace(4.5,2.175,70)
        x_inv_bt = []
        x_at = np.linspace(2.175,1.5,50)
        x_inv_at = []

        fit_to_plot_bt = []
        fit_to_plot_at = []
        for k in range(0,len(x_bt)):
            fit_to_plot_bt.append(Plotting_BCS_bt(x_bt[k],a0_fit,a1_fit,R0_fit,DeltaRs_fit))
            x_inv_bt.append(1/x_bt[k])
        for k in range(0,len(x_at)):
            fit_to_plot_at.append(Plotting_BCS_at(x_at[k],a0_fit,a1_fit,R0_fit,DeltaRs_fit))
            x_inv_at.append(1/x_at[k])

        step_plot_x = [x_inv_bt[-1],x_inv_at[0]]
        step_plot_y = [fit_to_plot_bt[-1], fit_to_plot_at[0]]

        #ax1.errorbar(Inv_Tdata_sep[i],Rs_sep[i], yerr=Rs_err_sep[i], marker='o', linestyle='none', markersize=4, color=colors[i], label=legend_entries[i])
        ax1.plot(Inv_Tdata_sep[i],Rs_sep[i], marker='o', linestyle='none', markersize=4, color=colors[i], label=legend_entries[i])
        ax1.step(step_plot_x, step_plot_y, marker='None', linestyle='--',color=colors[i])
        #ax1.plot(T_vals, Rs_for_T_vals[i], marker='None', linestyle='--',color=colors[i])
        #ax1.plot(Inv_Tdata_sep_bt, np.exp(result.best_fit[:len(Inv_Tdata_sep_bt)]), marker='None', linestyle='--',color=colors[i])
        #ax1.plot(Inv_Tdata_sep_at, np.exp(result.best_fit[len(Inv_Tdata_sep_bt):]), marker='None', linestyle='--',color=colors[i])
        ax1.plot(x_inv_bt, fit_to_plot_bt, marker='None', linestyle='--',color=colors[i])
        ax1.plot(x_inv_at, fit_to_plot_at, marker='None', linestyle='--',color=colors[i])
        #print(fit_to_plot_bt)
        #ax1.plot(Inv_Tdata_sep_at, Plotting_BCS_at(Tdata_sep_at,a0_fit,a1_fit,R0_fit,DeltaRs_fit), marker='None', linestyle='--',color=colors[i])
        #print(a0_fit, a1_fit, R0_fit, DeltaRs_fit)
        #print(Inv_Tdata_sep[i], Inv_T_vals)
        ax2.plot(Inv_Tdata_sep[i],((result.residual))*100/Rs_sep_ln[i], marker='o', markersize=3, color=colors[i], label=legend_entries[i])

    elif i < (len(colors)+len(shapes)):
        ax1.plot(Inv_Tdata_sep[i],Rs_sep[i], marker=shapes[i-len(colors)], linestyle='none', markersize=4, color='black', label=legend_entries[i])
        ax1.plot(Inv_Tdata_sep[i], result.best_fit, marker='None', linestyle='--',color='black')
        ax2.plot(Inv_Tdata_sep[i],((result.residual))*100/Rs_sep[i], marker=shapes[i-len(colors)], markersize=3, color='black', label=legend_entries[i])
    else:
        ax1.plot(Inv_Tdata_sep[i],Rs_sep[i], marker='o', linestyle='none', markersize=4, color='black', label=legend_entries[i])
        ax1.plot(Inv_Tdata_sep[i], result.best_fit, marker='None', linestyle='--',color='black')
        ax2.plot(Inv_Tdata_sep[i],((result.residual))*100/Rs_sep[i], marker='o', markersize=3, color='black', label=legend_entries[i])

file_fit_params.close()

for j in range(0,len(fixed_temps)):

    for k in range(0,len(FieldValues)):
        #print(FieldValues[k])
        Rss_val = Rev_Rs_correction(Rs_fixed_temps[j],FieldValues[k])

        Rss_fixed_temps[j].append(Rss_val)
        Q_fixed_temps[j].append((G_vals[cavity]*1e9)/Rss_val)
"""
for i in range(0, len(fixed_temps)):

    params = fixed_T_model.make_params(alpha=1, beta=1, gamma=1)
    result = fixed_T_model.fit(Rs_fixed_temps[i], params, B=FieldValues)

    str_beta_sign = '-'
    if result.best_values['beta'] >= 0:
        str_beta_sign = '+'

    str_gamma_sign = '-'
    if result.best_values['gamma'] >= 0:
        str_gamma_sign = '+'

    eq_str = ' Equation of fit line: ' + str(round(result.best_values['alpha'],2)) + ' ' +  str_beta_sign + ' ' +\
    str(abs(round(result.best_values['beta'],3))) + 'B' + ' ' + str_gamma_sign + ' ' + \
    str(abs(round(result.best_values['gamma'],5))) + r'B$^2$'

    ax3.plot(FieldValues, np.exp(Rs_fixed_temps[i]), marker='o', linestyle='none', markersize=4, label=(temps_legend[i]+eq_str), color=colors[i])
    ax3.plot(FieldValues, np.exp(result.best_fit), marker='None', linestyle='--',color=colors[i])
    #ax3.plot(FieldValues, Rs_for_T_vals[i], marker='None', linestyle='--',color=colors[i])
    ax4.plot(FieldValues, Q_fixed_temps[i], marker='o', linestyle='none', markersize=4, label=(temps_legend[i]), color=colors[i])
"""
for i in range(0, len(fixed_temps)):

    max_guess = Rs_fixed_temps[i][0]
    params = fixed_T_quad_model.make_params(R0q=20, gammaq=9)
    params['R0q'].set(min=0,max=(max_guess+0.5))
    params['gammaq'].set(min=0,max=15)
    result = fixed_T_quad_model.fit(Rs_fixed_temps[i], params, B=FieldValues, method='differential_evoluti', max_nfev=200000)

    str_gammaq_sign = '-'
    if result.best_values['gammaq'] >= 0:
        str_gammaq_sign = '+'

    eq_str = ' Equation of fit line: ' + str(round(result.best_values['R0q'],3)) + ' ' + str_gammaq_sign + ' ' +\
    str(round(result.best_values['R0q'],3)) + '*' + str(round(result.best_values['gammaq'],3)) + r'$(B/B_0)^2$'

    ax3.plot(FieldValues, (Rs_fixed_temps[i]), marker='o', linestyle='none', markersize=4, label=(temps_legend[i] + eq_str), color=colors[i])
    ax3.plot(FieldValues, (result.best_fit), marker='None', linestyle='--',color=colors[i])


#print("Rs T vals")
#print(Rs_fixed_temps)
#print()
#print(Rss_fixed_temps)
#format the x axis
x_formatter = FixedFormatter([r'10$^{-1}$', r'9$^{-1}$', r'8$^{-1}$', r'7$^{-1}$', r'6$^{-1}$', r'5$^{-1}$',
 r'4.5$^{-1}$', r'4.0$^{-1}$', r'3.5$^{-1}$', r'3.0$^{-1}$', r'2.5$^{-1}$', r'2.2$^{-1}$', r'2.0$^{-1}$', r'1.8$^{-1}$', r'1.7$^{-1}$'])
x_locator = FixedLocator([0.1, 1/9, 0.125, 1/7, 1/6, 0.2, 1/4.5, 0.25, 1/3.5, 1/3, 0.4, 1/2.2, 0.5, 1/1.8, 1/1.7])

ax1.set_yscale('log')
ax1.set_xlabel(r'Inverse Temperature [K$^{-1}$]',fontsize=14)
ax1.set_ylabel(r'R$_s[n\Omega]$',fontsize=14)
ax1.xaxis.set_major_formatter(x_formatter)
ax1.xaxis.set_major_locator(x_locator)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.vlines(x=0.45977, ymin=0, ymax=200, color='gray', linestyle='dashed', label = 'lambda point')
ax1.legend(title='Field Amplitude',fontsize=14)
ax1.set_title("Surface Resistance vs Inverse Temperature",fontsize=18)
ax1.grid(True)

ax2.set_xlabel(r'Inverse Temparature [K$^{-1}$]')
ax2.set_ylabel(r"Residual/R$_s$ [%]")
ax2.xaxis.set_major_formatter(x_formatter)
ax2.xaxis.set_major_locator(x_locator)
ax2.set_title(r"Percent Difference between observed R$_s$ and values calculated from fit parameters")
ax2.legend(title='Field Amplitude')
ax2.grid(True)

#ax3.set_yscale('log')
ax3.set_xlabel('RF field [mT]',fontsize=14)
ax3.set_ylabel(r'R$_s[n\Omega]$',fontsize=14)
ax3.legend(title='Temperature [K]',fontsize=14)
ax3.tick_params(axis='both', which='major', labelsize=14)
ax3.set_title('Total Surface Resistance vs RF Field for Fixed Temperatures')
ax3.grid(True)

ax4.set_xlabel('RF field [mT]',fontsize=14)
ax4.set_ylabel('Q',fontsize=14)
ax4.legend(title='Temperature [K]',fontsize=14)
ax4.tick_params(axis='both', which='major', labelsize=14)
#ax4.set_title('Total Surface Resistance vs RF Field for Fixed Temperatures')
ax4.grid(True)

plt.show()
