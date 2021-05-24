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

print("Some message about what this script is")
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

Qcol_str = '11'
Tcol_str = '24'
Eacccol_str = '9'
G_str = '120'
skiplines_str = '4'

skiplines = 0 #The default number of lines to skip at the beginning of the table is zero

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

for i in range(0,len(Eaccdata)):

    if i == 0:
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])
        continue

    if Eaccdata[i-1] < Eaccdata[i]:
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])
    else:

        #Do the fit here
        coef = np.polyfit(ramp_Eacc,ramp_Rss,3) #3rd degree polynomial
        #print(coef)
        fit = np.poly1d(coef)
        #print(fit)

        #Get ready to start the next ramp up fit
        ramp_Eacc.clear()
        ramp_Rss.clear()
        ramp_Eacc.append(Eaccdata[i])
        ramp_Rss.append(Rssdata[i])

print(ramp_Eacc)

print(len(Tdata), len(Rssdata), len(Eaccdata))
#plt.plot(Tdata, Rsdata)
#plt.show()
