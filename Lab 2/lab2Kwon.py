#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 15:26:27 2022
Modified on Fri Mar 11 15:26:27 2022
@author: Your name

Description
------------
"""

# Part 1 - Histogram the data

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.stats as st

fig = plt.figure(1,figsize=(6,6))
fig.clf()
axes = [fig.add_subplot(321),\
        fig.add_subplot(322),\
        fig.add_subplot(323),\
        fig.add_subplot(324),\
        fig.add_subplot(325),\
        fig.add_subplot(326)]

# You can use axes[0], axes[1], ....  axes[5] to make the six histograms.

# Your code goes here

dataset = np.loadtxt('refractionData.txt', skiprows=3)
#print(dataset)

alpha10 = dataset[0]
alpha20 = dataset[1]
alpha30 = dataset[2]
alpha40 = dataset[3]
alpha50 = dataset[4]
alpha60 = dataset[5]

axes[0].hist(alpha10, range=(-10,50), bins=15)
axes[0].set_title('alpha = 10 deg.',fontsize=12)
axes[0].set_xlabel('beta (deg.)')
axes[0].set_xticks([-10,0,10,20,30,40,50])

axes[1].hist(alpha20, range=(-10,50), bins=15)
axes[1].set_title('alpha = 20 deg.',fontsize=12)
axes[1].set_xlabel('beta (deg.)')
axes[1].set_xticks([-10,0,10,20,30,40,50])

axes[2].hist(alpha30, range=(-10,50), bins=15)
axes[2].set_title('alpha = 30 deg.',fontsize=12)
axes[2].set_xlabel('beta (deg.)')
axes[2].set_xticks([-10,0,10,20,30,40,50])

axes[3].hist(alpha40, range=(-10,50), bins=15)
axes[3].set_title('alpha = 40 deg.',fontsize=12)
axes[3].set_xlabel('beta (deg.)')
axes[3].set_xticks([-10,0,10,20,30,40,50])


axes[4].hist(alpha50, range=(-10,50), bins=15)
axes[4].set_title('alpha = 50 deg.',fontsize=12)
axes[4].set_xlabel('beta (deg.)')
axes[4].set_xticks([-10,0,10,20,30,40,50])

axes[5].hist(alpha60, range=(-10,50), bins=15)
axes[5].set_title('alpha = 60 deg.',fontsize=12)
axes[5].set_xlabel('beta (deg.)')
axes[5].set_xticks([-10,0,10,20,30,40,50])

plt.tight_layout()


#%%

# Part 2 - Table of measurements

# Your code goes here
def torad(x):
    return x*(np.pi/180)
def tosine(x):
    return np.sin(x)
#alpha in radians------------
alpha = np.array([10,20,30,40,50,60])
ralpha = torad(alpha)

#sine of ralpha---------------
sinralpha = tosine(ralpha)

#beta in radians
rbeta = torad(dataset)
#mean value of rbeta---------------
rmeanbeta = np.zeros(len(alpha))
for i in range(len(alpha)):
    rmeanbeta[i] = np.mean(rbeta[i])
#stdev of rbeta
stdrmeanbeta = np.zeros(len(alpha))
for i in range(len(alpha)):
    stdrmeanbeta[i] = np.std(rbeta[i], ddof=1)/np.sqrt(len(rbeta[0])) #ddof=1 corrects so its the same value as excel/matlab/etc


#sine of rmeanbeta-----------
sinrmeanbeta = tosine(rmeanbeta)

#propogate unc of sinrmeanbeta -----------------
stdsinrmeanbeta = np.zeros(len(alpha))
for i in range(len(alpha)):
    stdsinrmeanbeta[i] = np.sqrt((np.cos(rmeanbeta[i])*stdrmeanbeta[i])**2)

#print table (x ax alpha and y ax values)
print("+"*74)
titles = 'Table of measurements (in Radians)'
print(f"||{titles:^70}||")
print("+"*74)
headers = np.array(['Alpha','Sin(Alpha)','Beta (mean)','Sin(Beta)', 'Unc. of Sin(Beta)'])
#alpha at the top
print(f"||{headers[0]:^7}||{headers[1]:^12}||{headers[2]:^13}||{headers[3]:^11}||{headers[4]:^19}||")
print('+'*74)
for i in range(6):
    print(f"||{ralpha[i]:^7.4f}||{sinralpha[i]:^12.4f}||{rmeanbeta[i]:^13.4f}||{sinrmeanbeta[i]:^11.4f}||{stdsinrmeanbeta[i]:^19.4f}||")
print('+'*74)




#%%

# Part 3 - Snells law plot and fit

fig = plt.figure(2,figsize=(6,6))
fig.clf()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

# You can use ax1 and ax2 for the Snell's law plot and the chi squared plot.

# Your code goes here

#sin b vs sin a
ax1.errorbar(sinralpha, sinrmeanbeta ,yerr=stdsinrmeanbeta, fmt="ok", markersize=4)

def linfit(sina, n):
    return sina/n # n is index of refraction
p0 = 1.0
params, covar = curve_fit(linfit, sinralpha, sinrmeanbeta, p0, sigma=stdsinrmeanbeta, absolute_sigma=True)
fitn = params[0]
errn = np.sqrt(covar[0][0])
print(f'The best fit value of the index of refraction is {fitn:.4f} +/- {errn:.4f} .')

xcurve = np.linspace(0.0,0.9,81)
fitsinrbeta = linfit(sinralpha,fitn)
curvesinrbeta = linfit(xcurve,fitn)

ax1.plot(xcurve,curvesinrbeta, color='b')

ax1.set_title("Snell's Law", fontsize=16)
ax1.set_xlabel('sin($\\alpha$)', fontsize=13)
ax1.set_ylabel('sin($\\beta$)', fontsize=13)
ax1.set_xticks([0.0,0.2,0.4,0.6,0.8])
ax1.set_yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.5,0.6])

chisqterms = np.zeros(len(alpha))
for i in range(len(alpha)):
    chisqterms[i] = ((sinrmeanbeta[i]-fitsinrbeta[i])**2)/(stdsinrmeanbeta[i]**2)
chisq = chisqterms.sum()
v = 5 #n=6, 1 fitting param
pvalue = st.chi2.sf(chisq,v)
print(f'The chi^2 for the best fit line is {chisq:.4f}. With {v} degrees of freedom, the P-value is {pvalue:.4f}.')
ax1.text(0.3,0.5,f'n = {fitn:.3f} $\\pm$ {errn:.3f}')

# Part 4 - Chi squared plot

# Your code goes here

#chi2 as a function of index of refraction
def chi2vn(n):
    nchisqterms = np.zeros(len(alpha))
    nfitsinrbeta = linfit(sinralpha,n)
    for i in range(len(alpha)):
        nchisqterms[i] = ((sinrmeanbeta[i]-nfitsinrbeta[i])**2)/(stdsinrmeanbeta[i]**2)
    nchisq = nchisqterms.sum()
    return nchisq
nvals = np.linspace(1.45,1.53,51)
chi2vals = np.zeros(51)
for i in range(51):
    chi2vals[i] = chi2vn(nvals[i])
chi2min = chi2vn(fitn)
chi2max = chi2vn(fitn+errn)
ax2.plot(nvals,chi2vals, color='black', )
ax2.set_xlim(1.45,1.53)
ax2.set_ylim(2.5,5.0)

#lines
ax2.vlines(fitn,2.5,4,'black','--',linewidth=1.5)
ax2.vlines(fitn+errn,2.5,4,'black','--',linewidth=1.5)
ax2.vlines(fitn-errn,2.5,4,'black','--',linewidth=1.5)
ax2.hlines(chi2min,1.45,1.53,'black',linewidth=1.5)
ax2.hlines(chi2max,1.45,1.53,'black',linewidth=1.5)
ax2.set_title('$\\chi^2$ vs index of refraction', fontsize=16)
ax2.set_xlabel('index of refraction', fontsize=13)
ax2.set_ylabel('$\\chi^2$', fontsize=13)
ax2.text(1.48,4.5,f'n = {fitn:.3f} $\\pm$ {errn:.3f}')


plt.tight_layout()
######################################################
