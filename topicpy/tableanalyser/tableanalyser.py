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
import pandas as pd
from matplotlib import pyplot as plt
from topicpy.tacos_plot import scatterdense
from scipy import stats as st
from scipy import stats
import os


def geneinfo(genename, df, nfiles,metric='fpkm'):
    """Extimate mean and var for a ENSG gene

    Keyword arguments:
    :param genename: ENSG
    :param genedata: list or array-like object across tissues

    :return:  gene dictionary
    """
    #gene = 123
    #genename = df['gene'][gene]
    print("name: %s"%genename)
    if 'fpkm' in metric:
        genedata = np.array([fpkm for fpkm in df.loc[df['gene']==genename].loc[:,df.keys()[1:]].values.reshape(nfiles,1) if (fpkm>1e-1)&(fpkm<1e5)])
    else:
        genedata = np.array([c for c in df.loc[df['gene']==genename].loc[:,df.keys()[1:]].values.reshape(nfiles,1)])
    try:
        genemean = np.nanmean(genedata)
        genevariance = np.nanvar(genedata)
        geneoccurrence = float(len(genedata.nonzero()[0]))/len(genedata)
        genedict = {
            'name':genename,
            'avg' : genemean,
            'var': genevariance,
            'data' : genedata,
            'occ' : geneoccurrence
        }
        try:
            q = mg.getgenes(genename[:15], 'name,symbol,refseq.rna,type_of_gene,bp')[0]
            print("Descr: %s"%q['name'])
            print("Symbol: %s"%q['symbol'])
            genedict['type']=q['type_of_gene']
        except:
            genedict['type']='unknown'
            pass
        print("mean: %10.2f"%genemean)
        print("var: %10.2f"%genevariance)
        print("occ: %10.2f"%(geneoccurrence))
    except:
        return {}
    return genedict

def genedistr(genedict, bins = 50, ax = None, density=False, label='', save=True, metric='fpkm', logx=True, logy=True):
    """
    Plot distriution across tissues
    """
    maxfpkm = np.max(genedict['data'])
    width = float(maxfpkm) / bins
    _range = (0 - 0.5 * width, maxfpkm + 0.5 * width)
    if ax == None:
        fig = plt.figure(figsize=(18, 15))
        ax = fig.subplots()
    else:
        fig=ax.get_figure()
    n, bin_edges, _ = ax.hist(genedict['data'], lw=1.5, density=density, histtype='step', range=_range, bins=bins, label=label)
    ax.set_title(genedict['name'], fontsize=35)
    ax.set_xlabel('%s'%metric, fontsize=35)
    ax.set_ylabel('#', fontsize=35)
    if logy:
        ax.set_yscale('log')
    if logx:
        ax.set_xscale('log')
    if save:
        fig.savefig("plot/genes/%s_distr.pdf"%(genedict['name']))

def geneplot(genedict, metric='fpkm'):
    """
    Plot FPKM across tissues
    """
    fig = plt.figure(figsize=(15, 5))
    plt.plot(genedict['data'], 'ob')
    plt.title(genedict['name'], fontsize=35)
    plt.xlabel("sample", fontsize=35)
    plt.ylabel("%s"%metric, fontsize=35)
    plt.yscale('log')
    plt.ylim(ymin=1e-4)
    plt.show()
    fig.savefig("plot/genes/%s_data.pdf"%(genedict['name']))

def genecoord(genedict, means, variances, metric='fpkm'):
    """
    plot gene position in gobal plot
    """
    fig = plt.figure(figsize=(18,8))
    plt.scatter(means, variances)
    plt.scatter([np.average(genedict['data'])],[np.var(genedict['data'])], marker='x', c='r', s=90, label=genedict['name'])
    plt.xlabel("$<%s>$"%metric, fontsize=35)
    plt.ylabel("$\sigma^2_{%s}$"%metric, fontsize=35)
    plt.yscale('log')
    plt.xlim(5e-5,np.power(10,np.log10(means.max())+1))
    plt.ylim((variances[variances.nonzero()].min()/10,np.power(10,np.log10(variances.max())+1)))
    plt.legend()
    plt.show()
    fig.savefig("plot/genes/%s_coord.png"%(genedict['name']))

def discretize(column, nquartiles=10):
    quantiles = np.quantile(column,np.arange(0,1,1./nquartiles)[1:])
    c = np.digitize(column, quantiles) + 1
    c[column==0]=0
    return c

def discretize_df_columns(df):
    qdf = pd.DataFrame(index=df.index.values)
    for i,sample in enumerate(df.columns.values):
        print(sample,i)
        column = df.loc[:,sample].values[0]
        s = pd.Series(data = discretize(column), index=df.index.values, dtype=int)
        s.name=sample
        qdf.insert(0,s.name,s.values.round(0))
    return qdf

try:
    df_symbols= pd.read_csv("https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_pub_ensembl_id&col=md_ensembl_id&col=md_eg_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit", index_col=[0], sep='\t')
except:
    df_symbols = pd.read_csv("https://raw.githubusercontent.com/fvalle1/topicpy/develop/topicpy/tableanalyser/gene_conversions.txt", index_col=[0], sep='\t')

def get_symbol(ensg):
    '''
    convert ensg to symbol
    
    :param str ensg:
    '''
    if ensg[:15] in df_symbols['Ensembl gene ID'].values:
        return df_symbols[df_symbols['Ensembl gene ID']==ensg[:15]]['Approved symbol'].values[0]
    else:
        return ''

def get_ensg(description):
    '''
    convert descr to ensg
    
    :param str description:
    '''
    if description in df_symbols['Approved symbol'].values:
        return df_symbols[df_symbols['Approved symbol']==description]['Ensembl gene ID']
    else:
        return ''

def plotvarmen(means, variances, ax = None, normalisation_str = "counts", poisson_limit=1, colorbar=False, **kwargs):
    x_lin = np.logspace(np.log10(np.nanmin(means)),np.log10(np.nanmax(means)), dtype=float,num=50)
    if ax is None:
        fig=plt.figure(figsize=(18,15))
        ax=fig.subplots()
    ax.plot(x_lin[x_lin>=poisson_limit-0.5],np.power(x_lin[x_lin>=poisson_limit-0.5],2)*1/poisson_limit, 'b--', lw=5, label="$m_g^2$ (Taylor)")
    ax.plot(x_lin[x_lin<=poisson_limit+0.5],x_lin[x_lin<=poisson_limit+0.5], 'r--', lw=5, label="$m_g$ (Poisson)")

    scatterdense(means, variances, ax=ax, label='data', colorbar=colorbar, c_title="Number of genes", **kwargs)
    ax.tick_params(labelsize=35, width=8, length=20)
    ax.set_xlabel("Mean expression level, $m_g$", fontsize=35)
    ax.set_ylabel("Variance of expression level, $v_g$", fontsize=35)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(np.nanmin(means[means.nonzero()])/5,np.power(10,np.log10(np.nanmax(means))+1))
    ax.set_ylim((np.nanmin(variances[variances.nonzero()])/10,np.power(10,np.log10(np.nanmax(variances))+1)))
    ax.legend(fontsize=35)

def plotcv2mean(means, variances, ax=None, normalisation_str = "counts", poisson_limit=1, colorbar=False, **kwargs):
    x_lin = np.logspace(np.log10(np.nanmin(means)),np.log10(np.nanmax(means)), dtype=float,num=50)
    cv2 = np.array([variances[i]/(np.power(mean,2)) for i,mean in enumerate(means) if mean>0])
    scatterdense(means[means>0], cv2,ax=ax, label="data", colorbar=colorbar, c_title="Number of genes", **kwargs)
    if ax is None:
        fig=plt.figure(figsize=(18,15))
        ax=fig.subplots()
    ax.plot(x_lin[x_lin>=poisson_limit-0.5],[1/poisson_limit for _ in x_lin[x_lin>=poisson_limit-0.5]], 'b--', lw=5, label="%.1f (Taylor)"%(1/poisson_limit))
    ax.plot(x_lin[x_lin<=poisson_limit+0.5],1./x_lin[x_lin<=poisson_limit+0.5], 'r--', lw=5, label="$m_g^{-1}$ (Poisson)")

    #plt.plot(x_lin, [nfiles-1 for _ in x_lin], color='cyan', ls='--', lw=3.5, label='bound')
    ax.tick_params(labelsize=35, width=8, length=20)
    ax.set_ylabel("Coefficient of variation squared, $CV^2_g$", fontsize=35)
    ax.set_xlabel("Mean expression level, $m_g$", fontsize=35)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(labelsize=20)
    ax.set_xlim(means[means>0].min()*0.9,np.power(10,np.log10(np.nanmax(means))+1))
    ax.set_ylim((cv2[cv2>0].min()/10,np.power(10,np.log10(np.nanmax(cv2)+1))))
    ax.legend(fontsize=24)


def plotoversigmacv2(means,variances, ax=None, normalisation_str = "counts", how_many_sigmas=3):
    x_lin = np.logspace(np.log10(means[means.nonzero()].min()),np.log10(means[means.nonzero()].max()), dtype=float,num=50)
    cv2 = np.array([variances[i]/(np.power(mean,2)) for i,mean in enumerate(means) if mean>0])
    if ax is None:
        fig=plt.figure(figsize=(18,15))
        ax=fig.subplots()
    ax.scatter(means[means.nonzero()], cv2, c='b')
    ax.tick_params(labelsize=35, width=8, length=20)
    ax.plot(x_lin[-30:],[1e-1 for _ in x_lin[-30:]], 'g-', lw=3.5, label='$1$ (Taylor)')
    ax.plot(x_lin[:30],1./x_lin[:30], 'r-', lw=3.5, label='$<%s>^{-1}$ (Poisson)'%normalisation_str)

    #plt.plot(x_lin, [nfiles-1 for _ in x_lin], color='cyan', ls='--', lw=3.5, label='bound')

    log_bins_for_x = np.logspace(np.log10(means[means.nonzero()].min()),6,25)

    bin_means, bin_edges, _ = st.binned_statistic(means[means.nonzero()], cv2, statistic='median', bins=log_bins_for_x)
    ax.scatter((bin_edges[:-1]+bin_edges[1:])/2.,bin_means, marker='x', color='orange', label='binned median')

    bin_sigmas,  _, _ = stats.binned_statistic(means[means.nonzero()], cv2, statistic=np.std, bins=log_bins_for_x)
    ax.plot((bin_edges[:-1] + bin_edges[1:])/2, bin_means+bin_sigmas*how_many_sigmas, lw=3, color='yellow', label='binned average + $%d\sigma$'%how_many_sigmas)


    ax.set_xlabel("$<%s>$"%normalisation_str, fontsize=35)
    ax.set_ylabel("$cv^2$", fontsize=35)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(means[means.nonzero()].min()/5,np.power(10,np.log10(means.max())+1))
    ax.set_ylim((cv2[cv2.nonzero()].min()/10,np.power(10,np.log10(cv2.max())+1)))
    ax.legend(fontsize=35)

def plotoverpoints(means, variances, over_plot, ax=None, normalisation_str = "counts", how_many_sigmas=3):
    if ax is None:
        fig=plt.figure()
        ax=fig.subplots()
    x_lin = np.logspace(np.log10(means[means.nonzero()].min()),np.log10(means[means.nonzero()].max()), dtype=float,num=50)
    cv2 = np.array([variances[i]/(np.power(mean,2)) for i,mean in enumerate(means) if mean>0])
    ax.scatter(means[means.nonzero()], cv2, c='b')
    ax.scatter(over_plot.T[0],over_plot.T[1], color='cyan')

    ax.plot(x_lin[-30:],[1e-1 for _ in x_lin[-30:]], 'g-', lw=3.5, label='(Taylor)')
    ax.plot(x_lin[:30],1./x_lin[:30], 'r-', lw=3.5, label='$<%s>^{-1}$ (Poisson)'%normalisation_str)
    #plt.plot(x_lin, [nfiles-1 for _ in x_lin], color='cyan', ls='--', lw=3.5, label='bound')

    log_bins_for_x = np.logspace(np.log10(means[means.nonzero()].min()),5,25)

    #bin_means, bin_edges, _ = st.binned_statistic(means[means.nonzero()], cv2, statistic='median', bins=log_bins_for_x)
    #plt.scatter((bin_edges[:-1]+bin_edges[1:])/2.,bin_means, marker='x', color='orange', label='binned median')
    bin_means, bin_edges, _ = st.binned_statistic(means[means.nonzero()], cv2, statistic='median', bins=log_bins_for_x)
    bin_sigmas,  _, _ = stats.binned_statistic(means[means.nonzero()], cv2, statistic=np.std, bins=log_bins_for_x)
    ax.hlines(bin_means+bin_sigmas*how_many_sigmas,bin_edges[1:], bin_edges[:-1], lw=3, color='yellow', label='binned average + $%d\sigma$'%how_many_sigmas)

    ax.tick_params(labelsize=35, width=8, length=20)
    ax.set_xlabel("$<%s>$"%normalisation_str, fontsize=35)
    ax.set_ylabel("$cv^2$", fontsize=35)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(means[means.nonzero()].min()/5,np.power(10,np.log10(means.max())+1))
    ax.set_ylim((cv2[cv2.nonzero()].min()/10,np.power(10,np.log10(cv2.max())+1)))
    ax.legend(fontsize=35)
    plt.show()

def getovergenes(df_mv, func, method='sampling', distance=10, how_many_sigmas=3, knee=100):
    '''
    :param method: can be 'sampling', 'sigma'
    '''
    over = []
    if 'sigma' in method:
        log_bins_for_x = np.logspace(np.log10(df_mv['mean'][df_mv['mean']!=0].min()),5,25)
        cv2=df_mv['variance']/df_mv['mean']/df_mv['mean']
        bin_means,bin_edges ,_ = st.binned_statistic(df_mv['mean'], cv2, statistic='mean', bins=log_bins_for_x)
        bin_sigmas, _ ,_ = st.binned_statistic(df_mv['mean'], cv2,statistic='std', bins=log_bins_for_x)
    for g in df_mv.index:
        subdf = df_mv.loc[g,:]
        mean = float(subdf['mean'])
        var = subdf['variance']
        if mean> 1e5 or mean< 1e-3:
            continue
        if 'sampling' in method:
            r = func(mean,float(var)/mean/mean, knee=knee, distance=distance)
        elif 'sigma' in method:
            r = func(mean,float(var)/mean/mean, bin_means=bin_means, bin_sigmas=bin_sigmas, bin_edges=bin_edges,how_many_sigmas=how_many_sigmas)
        if r[4]:
            over.append(g)
    print("found %d highly variable genes"%len(over))
    over_plot = []
    for g in over:
        subdf = df_mv.loc[g,]
        mean = subdf['mean']
        var = subdf['variance']
        occ = subdf['occurrence']
        cv2 = float(var)/mean/mean
        over_plot.append([mean,cv2,occ])
    over_plot=np.array(over_plot)
    return (over,over_plot)

def get_mean_cv2_sampling(mean, cv2, knee=1., distance=10):
    if mean < knee:
        return(mean, cv2, -1, -1, cv2 > distance+1./mean)
    else:
        return(mean, cv2, -1, -1, cv2 > 1e-1+distance)


def get_mean_cv2_sigma(mean, cv2, how_many_sigmas=3):
    bin_i = 0
    for i in range(len(bin_edges[:-1])):
        if mean<=bin_edges[i+1] and mean > bin_edges[i]:
            bin_i = i
            break
    return(mean, cv2, bin_means[bin_i], bin_sigmas[bin_i], cv2>(bin_means[bin_i]+how_many_sigmas*bin_sigmas[bin_i]))

def scalinglawsandoverexpressed(working_dir, normalisation_str = "counts", method='sampling', how_many_sigmas=3, distance=10):
    os.chdir(working_dir)
    df = pd.read_csv(("mainTable.csv"), index_col=[0])
    print(df.info())
    ngenes = len(df.index)
    nfiles = len(df.columns)
    print("genes:%d\trealizations:%d"%(ngenes,nfiles))
    df_mv = pd.read_csv("meanVariances.csv", index_col = [0])
    df_mv.dropna(axis=0,how='any',inplace=True)
    print(df_mv.info())
    means = df_mv['mean'].values
    variances = df_mv['variance'].values
    occurrences = np.array(df_mv['occurrence'].values, dtype=float)
    abundances = pd.read_csv("A.dat", header=None).values
    #vm
    fig=plt.figure(figsize=(15,8))
    ax=fig.subplots()
    plotvarmen(means, variances, ax=ax, normalisation_str=normalisation_str)
    fig.savefig("varmean_loglog.png")
    #cvm
    fig=plt.figure(figsize=(15,8))
    ax=fig.subplots()
    plotcv2mean(means, variances, ax=ax, normalisation_str=normalisation_str)
    fig.savefig("cvmean_loglog.png")
    #cvover
    fig=plt.figure(figsize=(15,8))
    ax=fig.subplots()
    plotoversigmacv2(means, variances, ax=ax, normalisation_str=normalisation_str, how_many_sigmas=how_many_sigmas)
    fig.savefig("cvmean_loglog_%dsigma.png"%how_many_sigmas)
    if 'sampling' in method:
        over, over_plot = getovergenes(df_mv,get_mean_cv2_sampling, method='sampling', distance=distance)
    elif 'sigma' in method:
        over, over_plot = getovergenes(df_mv,get_mean_cv2_sigma, method='sigma')
    #cvoverpoints
    fig=plt.figure(figsize=(15,8))
    ax = fig.subplots()
    plotoverpoints(means, variances, over_plot, ax=ax, how_many_sigmas=how_many_sigmas)
    fig.savefig("cvmean_loglog_oversigma.png")
    df[df.index.isin(over)].dropna().to_csv("mainTable_over.csv",index=True, header=True)
