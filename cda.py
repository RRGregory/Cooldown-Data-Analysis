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
import matplotlib.pyplot as plt
import numpy as np

print("Running the script")
print()

#Get file path and name from user
#file_path = input("Enter the full path and name of your input data file: ")
file_path = '/Users/ruthgregory/Documents/SRF/Data/before_bake/QWR_217MHz_cooldown_10uT_2019_12_06_08h30_QoData.txt'
file = open(file_path, "r")

#Get information about the input data table from user
#Qcol_str = input("Enter the column number of your Quality Factor data: ")
#Tcol_str = input("Enter the column number of your temperature data: ")
#Eacccol_str = input("Enter the column number of your accelerating field data: ")
#G_str = input("Enter the Geometric factor: ")
#skiplines_str = input("Optional input: Enter the number of lines to skip at the beginning of the table, or to not skip any lines, press enter: ")

cavity = 2 #Value for gemetric factor, freqeuncy, and determining which betas to use
Qcol_str = '11'
Tcol_str = '24'
Eacccol_str = '9'
G_str = '120'
skiplines_str = '4'

FieldValues = [10,20,30,40,50,60,70]# #field amplitude levels, array input parameter

skiplines = 0 #The default number of lines to skip at the beginning of the table is zero

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

#Convert the user input values to ints
if skiplines_str != "":
    skiplines = int(skiplines_str)

Qcol = int(Qcol_str)
Tcol = int(Tcol_str)
Eacccol = int(Eacccol_str)
G = int(G_str)

#Subtract 1 because Python indexing starts at 0
Qcol -= 1
Tcol -= 1
Eacccol -= 1

#Create empty lists to hold the data
Qdata = []
Tdata = []
Eaccdata = []
Rssdata = [] #This will be the non-corrected Rs* data: G/Qo

def Rs_correction(Bp, Rparam):
    return beta3_all[cavity]*Rparam[0]*Bp**3 + beta2_all[cavity]*Rparam[1]*Bp**2 + beta1_all[cavity]*Rparam[2]*Bp + beta0_all[cavity]*Rparam[3]


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
        Qdata.append(float(columns[Qcol]))
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

#Get the non-corrected Rs* data
for Q in Qdata:
    Rssdata.append(G/Q)


ramp_Eacc = []
ramp_Rss = []
Rs = [] #Array for corrected Rs data

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
        fit = np.poly1d(coef)
        #print(fit)

        for j in range(0,len(ramp_Rss)):
            Rs.append(Rs_correction(ramp_Eacc[j], coef))
        print("last ramp", ramp_Eacc[0])

    if Eaccdata[i-1] < Eaccdata[i]:
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])
    else:

        #Do the fit here
        coef = np.polyfit(ramp_Eacc,ramp_Rss,3) #3rd degree polynomial
        #print(coef)
        fit = np.poly1d(coef)
        #print(fit)

        for j in range(0,len(ramp_Rss)):
            Rs.append(Rs_correction(ramp_Eacc[j], coef))
        print("some ramp", ramp_Eacc[0])

        #Get ready to start the next ramp up fit
        ramp_Eacc.clear()
        ramp_Rss.clear()
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])

#print(ramp_Eacc)
print(len(Tdata), len(Rssdata), len(Eaccdata), len(Rs))


#Make list of file names to split the data by field amplitude
file_names = []
for i in range(0,len(FieldValues)):

    name = 'temp_Rs_' + str(FieldValues[i]) + 'uT.txt'
    file_names.append(name)
    #print(file_names[i])

for i in range(0,len(Eaccdata)):

    for j in range(0,len(FieldValues)):

        if (Eaccdata[i] < FieldValues[j]+0.5) and (Eaccdata[i] > FieldValues[j]-0.5):

            entry = str(Tdata[i]) + "," + str(Rs[i]) + "\n"
            file = open(file_names[j], "a")
            file.write(entry)
            file.close()
