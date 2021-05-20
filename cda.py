"""
 *========================================================================================
 *
 *       File name:  cda.py
 *
 *    Description:  Takes Q vs T data collected during a cooldown and extracts
 *                  Rbcs and Rres as a function of field amplitude
 *
 *        Version:  1.0
 *        Created:  19/05/2021
 *       Revision:  none
 *       Compiler:  Python 3.8
 *
 *         Author:  Ruth Gregory
 *          Email:  ruthannrgregory@gmail.com
 *   Organization:  TRIUMF
 *
 * =======================================================================================
 """
#Comments from Philipp's script:

# input parameters shoudld be:
# - data file with T,field amplitude,  and Q0 as columns
# - which columns the T, Eacc, Q0 data is in
# - geometric factor G
# - RF frequency (watch unit, MHz in this code)
# - beta factors for field distribution correction (array of four floating point numbers)
# - measured field amplitude
# - SWR for uncertainty
# - output file(s) that contain:
    # - fit results and uncertainties
    # - fit quality parameter (R^2..)

file_name = "QWR_217MHz_cooldown_10uT_2019_12_06_08h30_QoData.txt"
file_path = "/Users/ruthgregory/Documents/SRF/Data/before_bake"

file = open("/Users/ruthgregory/Documents/SRF/Data/before_bake/QWR_217MHz_cooldown_10uT_2019_12_06_08h30_QoData.txt", "r")
#data = file.read()
Qdata = []
Rsdata = []
G = 120

line_cnt = 1 #line count, it will start with line 1
for line in file:
    line = line.strip()
    columns = line.split()
    if columns == []:
        line_cnt += 1
        continue

    try:
        print(columns[10])
        Qdata.append(float(columns[10]))
    except IndexError:
        print("Warning: index error on line ", line_cnt, "in Q column, skipping line", line_cnt)
    except ValueError:
        line_cnt+=1
        continue

    line_cnt+=1

file.close()

for Q in Qdata:
    Rsdata.append(G/Q)

print(Rsdata)
