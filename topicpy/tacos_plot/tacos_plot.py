#  Copyright (c) 2020 fvalle
#
#  Permission is hereby granted, free of charge, to any person
#  obtaining a copy of this software and associated documentation
#  files (the "Software"), to deal in the Software without
#  restriction, including without limitation the rights to use,
#  copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following
#  conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#  OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
from scipy.interpolate import interpn
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize



def scatterdense(x, y, ax=None, nbins=80, colorbar=False, c_title="density", **kwargs) -> None:
    """
    Plot a scatter plot with density
    
    :param x: list of first variable
    :param y: list of second variable
    :param ax: matplotlib.pyplot.Axes add here or create a new one
    :param nbins: number of bins
    :param colorbar: plot colorbar
    :param c_title: color bar title
    :param **kargs: arguments to be passed to matplotlib.pyplot.scatter
    """
    xmin = np.log10(min(x[x>0]))
    xmax = np.log10(max(x))
    ymin = np.log10(min(y[y>0]))
    ymax = np.log10(max(y))

    xbins = np.logspace(xmin, xmax, nbins) # <- make a range from 10**xmin to 10**xmax
    ybins = np.logspace(ymin, ymax, nbins) # <- make a range from 10**ymin to 10**ymax
    data , x_e, y_e = np.histogram2d(x, y, bins = (xbins, ybins))
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    if ax is None:
        fig=plt.figure()
        ax=fig.subplots()
    ax.scatter(x, y, c=z, **kwargs)
    if colorbar:
        cbar = ax.get_figure().colorbar(cm.ScalarMappable(norm=Normalize(vmin=1, vmax=max(z)), cmap="viridis"), ax=ax)
        cbar.ax.set_ylabel(c_title, fontsize=35)
