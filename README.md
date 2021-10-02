Comments from Philipp's script:

 input parameters shoudld be:
 - data file with T,field amplitude,  and Q0 as columns
 - which columns the T, Eacc, Q0 data is in
 - geometric factor G
 - RF frequency (watch unit, MHz in this code)
 - beta factors for field distribution correction (array of four floating point numbers)
 - measured field amplitude
 - SWR for uncertainty
 - output file(s) that contain:
     - fit results and uncertainties
     - fit quality parameter (R^2..)
------------------------------------------------------------------------------------------------------------
## Input Parameters

This script makes Python lists containing Temperature, field amplitude, and quality factor (Q<sub>0</sub>) data.
Note: the data input methodology is not really developed at this point, this will be a last step.


## Converting Rs* to Rs Data

At the beginning of this code there are lists of beta factors, described [here](https://journals.aps.org/prab/abstract/10.1103/PhysRevAccelBeams.21.122001). The beta factors used here have been calculated for six different RF modes: two for a QWR and four for a HWR. These beta factors will be used for converting the Rs* data to Rs data.

To start, this code calculates the non-corrected Rs values (denoted Rs*) by taking the quotient of the geometry factor (G) and the quality factor (Q<sub>0</sub>) for every data point. These values are stored in a Pyhton list object which has the same indexing as the lists containing the Temperature, Eacc, and Q<sub>0</sub> data.

The cooldown data can be separated into "ramp ups" of the cavity field. When taking data the cavity field will typically be increased from 10 to 60-100 mT, with one data point collected for every 10 mT. For each of these ramp ups a 3rd degree polynomial fit is made from the field dependent Rs* values and their corresponding peak fields (Rs*(B<sub>p</sub>). These values are fit to a polynomial function of the form:

<img src="https://render.githubusercontent.com/render/math?math=y_{uncorrected} = ax^3%20%2B%20bx^2%20%2B%20cx%20%2B%20d">

The coefficients of this polynomial are extracted and stored in an array. Finally, to calculate the field corrected Rs values for a given ramp up, the B<sub>p</sub> values from said ramp up are passed into a function exactly like the fit function calculated earlier, but with the polynomial coefficients multiplied by their appropriate beta factors:

<img src="https://render.githubusercontent.com/render/math?math=y_{corrected} = a\beta_3x^3%20%2B%20b\beta_2x^2%20%2B%20c\beta_1x%20%2B%20d\beta_0">

These corrected Rs values are then stored in a list with the same indexing as the other lists.

## Uncertainties

The uncertainties in the peak field (Bp) are:  <img src="https://render.githubusercontent.com/render/math?math=\deltaB_p = \frac{B_p(SWR-1)}{4}">

The uncertainties in the non-corrected surface resistance (Rs*) are:  <img src="https://render.githubusercontent.com/render/math?math=\deltaR_s^* = \frac{R_s^*(SWR-1)}{2}">

The uncertainties in the field corrected Rs values are given by 
![image](https://user-images.githubusercontent.com/19824886/135696597-99a407ae-4eef-4388-9319-d93a5f42ed62.png)

where in this code the  <img src="https://render.githubusercontent.com/render/math?math=\alpha_i"> values are just 0,1,2,3. 

The fitting described above is done with a least squares polynomial fit. The uncertainty in the fit parameters(  <img src="https://render.githubusercontent.com/render/math?math=r_{\alphai}"> ) was obtained from the diagonal of the covariance matrix as described in this paper: https://ipnpr.jpl.nasa.gov/progress_report/42-122/122E.pdf 

## Other notes and possible future improvements
- An unresolved issue is that sometimes the input data does not have all strictly increasing field amplitude values for each ramp up. For example the experimenters might be trying to take a measurement at 30 MV/m, overshoot and go to 33 MV/m, and then go back and take a data point at 30 MV/m. Then the data column for accelerating field data would have the series 10,20,33,30,40,50,60,70. Since the script looks for increasing values (within 0.5 mT) when it splits the data into ramp-ups, these would be split as [10,20,33] and [30,40,50,60,70], leading to less than ideal polyfits when doing the Rs* to Rs conversions described above. If you see the warning "sys:1: RankWarning: Polyfit may be poorly conditioned" when you run the script, this is probably what happened.
