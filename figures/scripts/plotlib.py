#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import statsmodels.api as sm
import math

def identity_line(ax=None, ls='--', *args, **kwargs):
    '''Draws identity line'''
    
    ax = ax or plt.gca()
    identity, = ax.plot([], [], ls=ls, *args, **kwargs)
    def callback(axes):
        low_x, high_x = ax.get_xlim()
        low_y, high_y = ax.get_ylim()
        low = min(low_x, low_y)
        high = max(high_x, high_y)
        identity.set_data([low, high], [low, high])
    callback(ax)
    ax.callbacks.connect('xlim_changed', callback)
    ax.callbacks.connect('ylim_changed', callback)
    return ax

def plot_line(intercept, slope, xmin, xmax, c, ax, ls='solid'):
    x = np.array([xmin, xmax])
    y = intercept + (slope * x)
    ax.plot(x, y, c=c, ls=ls)
    
def make_xtick_lab(xticks, xmaj):
    xlab = [''] * len(xticks)
    for i, x in enumerate(xticks):
        if x in xmaj:
            if (x > 1000) or (x < 0.001):
                factor, pwr = format(x, '0.0e').split('e')
                factor = int(factor); pwr = int(pwr)
                #pwr = int(fmt_str[1])
                #factor = int(math.floor(x / 10**pwr))
                if(factor == 1):
                    xlab[i] = r'$10^{:d}$'.format(pwr)
                else:
                    xlab[i] = r'${:d} \times 10^{:d}$'.format(factor, pwr)
            elif x < 1:
                xlab[i] = format(x, '0.1f')
            else:
                xlab[i] = format(x, '0.0f')
    return xlab

def make_pgls_legend(lam=None, pval=None, r2=None):
    outstr = ''
    if lam is not None:
        outstr += r'$\lambda = {:0.2f}$'.format(lam) + '\n'
    if pval is not None:
        if pval > 1e-3:
            outstr += r'$p = {:0.3f}$'.format(pval) + '\n'
        else:
            # Use scientific notation
            x, pwr = '{:0.1e}'.format(pval).split('e')
            pwr = '{:d}'.format(int(pwr)) 
            
            outstr += r'$p = {:s} \times 10^{{{}}}$'.format(x, pwr) + '\n'
    if r2 is not None:
        outstr += r'$r^2 = {:0.2f}$'.format(r2)
    return outstr

def get_CI(fit, x, y, bootstrap=False):
    xmin = min(x.to_numpy()[:,1])
    xmax = max(x.to_numpy()[:,1])
    xPred = np.column_stack([np.repeat(1.0, 100), np.linspace(xmin, xmax, 100)])
    if bootstrap:
        nBoot = 500
        n = y.shape[0]
        bootSamp = np.random.randint(0, n, size=(nBoot*n)).reshape((n, nBoot))
        bootCoef = np.array([sm.OLS(y.iloc[bootSamp[:,i]], x.iloc[bootSamp[:,i]]).fit().params for i in range(nBoot)])
        
        ci = np.quantile(np.array([np.sum(np.multiply(bootCoef[i,:], xPred), axis=1) for i in range(nBoot)]), q=(0.025, 0.975), axis=0)
        ci = ci.transpose()
    else:
        ci = fit.get_prediction(xPred).conf_int(alpha=0.05)

    outDict = {'xPred':xPred, 'ci':ci}
    
    return outDict

def squarify(n):
    nr = round(math.sqrt(n))
    nc = math.ceil(n / nr)
    return (nr, nc)

def make_log_xticks(xmin, xmax):
    
    for i in range(10):
        xstart = (i * 10**math.floor(xmin))
        if xstart >= 10**xmin:
            break

    for i in range(10):
        xend = (i * 10**math.floor(xmax))
        if xend > 10**xmax:
            xend = ((i-1) * 10**math.floor(xmax))
            break

    tick = xstart + 0
    interval = 10 ** math.floor(math.log10(xstart))
    xticks, xmaj = [], []
    while tick <= xend:
        xticks.append(tick)
        
        if (tick == 1) or (tick == xstart) or (math.log10(tick) % 1 == 0) or ((tick+interval > xend) and ((tick - xmaj[-1]) / (interval*10) > 0.25)):
            xmaj.append(tick)

        tick += interval
        if(tick == interval*10):
            interval *= 10

    # If xticks is too long, just mark off powers of 10
    if(len(xticks) > 20):
        for i in xticks[1:-1]:
            if (math.log10(i) % 1) != 0:
                xticks.remove(i)
        xmaj = xticks[:]

        if ((xmaj[-1] - xmaj[-2]) / (10**(math.floor(math.log10(xmaj[-1]))+1)) < 0.25):
            rm_x = xmaj[-1]
            xmaj.remove(rm_x)
            xticks.remove(rm_x)

    return (xticks, xmaj)

# if __name__ == '__main__':
#     
#     
