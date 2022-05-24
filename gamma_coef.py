import matplotlib.pyplot as plt
import matplotlib
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

Inv_T = [1/2.1,1/2.5,1/3,1/3.5,1/4]
Inv_T100 = [1/2.1,1/2.5,1/3,1/3.5]

Rss_10mT = [3.63693E-09,4.73795E-09,9.09831E-09,1.36173E-08,2.19047E-08]
Rss_20mT = [3.56424E-09,4.84463E-09,9.85311E-09,1.50207E-08,2.37802E-08]
Rss_30mT = [3.67529E-09,5.03585E-09,1.11244E-08,1.71066E-08,2.75039E-08]
Rss_40mT = [3.82534E-09,5.48774E-09,1.27455E-08,2.0341E-08,3.29903E-08]
Rss_50mT = [4.12576E-09,5.98573E-09,1.50816E-08,2.45933E-08,4.01611E-08]
Rss_60mT = [4.51305E-09,6.74646E-09,1.81328E-08,3.02269E-08,4.98458E-08]
Rss_70mT = [5.04564E-09,7.78411E-09,2.21434E-08,3.75767E-08,6.26232E-08]
Rss_80mT = [6.07087E-09,9.47229E-09,2.75771E-08,4.72047E-08,7.88133E-08]
Rss_90mT = [6.83637E-09,1.09874E-08,3.36447E-08,5.6314E-08,9.58062E-08]
Rss_100mT = [7.68135E-09,1.42754E-08,3.90742E-08,6.91596E-08]

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)

fig1, ax1 = plt.subplots(nrows=1, ncols=1)

x_formatter = FixedFormatter([r'10$^{-1}$',r'4.0$^{-1}$', r'3.5$^{-1}$', r'3.0$^{-1}$', r'2.5$^{-1}$', r'2.1$^{-1}$'])
x_locator = FixedLocator([0.1,0.25, 1/3.5, 1/3.0, 1/2.5, 1/2.1])

ax1.plot(Inv_T, Rss_10mT, marker='o', linestyle='none', markersize=6, label='10')
ax1.plot(Inv_T, Rss_20mT,marker='o', linestyle='none', markersize=6, label='20')
ax1.plot(Inv_T, Rss_30mT,marker='o', linestyle='none', markersize=6, label='30')
ax1.plot(Inv_T, Rss_40mT,marker='o', linestyle='none', markersize=6, label='40')
ax1.plot(Inv_T, Rss_50mT,marker='o', linestyle='none', markersize=6, label='50')
ax1.plot(Inv_T, Rss_60mT,marker='o', linestyle='none', markersize=6, label='60')
ax1.plot(Inv_T, Rss_70mT,marker='o', linestyle='none', markersize=6, label='70')
ax1.plot(Inv_T, Rss_80mT,marker='o', linestyle='none', markersize=6, label='80')
ax1.plot(Inv_T, Rss_90mT, marker='o', linestyle='none', markersize=6,label='90')
ax1.plot(Inv_T100, Rss_100mT,marker='o', linestyle='none', markersize=6, label='100')
ax1.legend(title='Peak Field [mT]',fontsize=20)
ax1.set_xlabel(r'Inverse Temperature [K$^{-1}$]',fontsize=24)
ax1.set_ylabel(r'R$_s*[\Omega]$',fontsize=24)
ax1.tick_params(axis='both', which='major', labelsize=24)
ax1.xaxis.set_major_formatter(x_formatter)
ax1.xaxis.set_major_locator(x_locator)
ax1.grid(True)

plt.show()
