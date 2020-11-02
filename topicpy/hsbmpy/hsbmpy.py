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
import os

import pandas as pd
import numpy as np
from sklearn.utils import shuffle
from matplotlib import pyplot as plt
import sys
import seaborn as sns

sns.set()
sns.set_context("paper")
from sklearn import metrics

# get colors from https://medialab.github.io/iwanthue/ or artenatevly from http://phrogz.net/css/distinct-colors.html
colors_cycle = ["#a257d4",
                "#e090bf",
                "#64c9a3",
                "#4b68ae",
                "#dc8c2f",
                "#cd41a7",
                "#d9344f",
                "#bc599a",
                "#afa1e8",
                "#48c1d8",
                "#b54545",
                "#919233",
                "#9a78be",
                "#59602a",
                "#4e8e2c",
                "#9db935",
                "#9b563c",
                "#e482df",
                "#5995d3",
                "#6a5198",
                "#b05f84",
                "#b563c3",
                "#5f6b18",
                "#a55c21",
                "#5754c2",
                "#277257",
                "#4f9b5e",
                "#8b6b29",
                "#b8381c",
                "#ad2f62",
                "#97ba6d",
                "#45c37c",
                "#5fc250",
                "#8c4c7b",
                "#e06e87",
                "#e2672a",
                "#db7756",
                "#974858",
                "#35743b",
                "#bbaf6c",
                "#8c4099",
                "#e44586",
                "#ed5c4c",
                "#389c84",
                "#cfae3d",
                "#eda377",
                "#778749",
                "#c5935a",
                "#de8784",
                "#757eec"]

def get_next_color():
    for color in colors_cycle:
        yield color

color_iterator = get_next_color()

def plot_cluster_composition(fraction_sites, directory, level, normalise=False, label='primary_site', shuffled=False,
                             algorithm='topsbm'):
    sns.set(font_scale=0.8)
    df_clusters = pd.read_csv("%s/%s/%s_level_%d_clusters.csv" % (directory, algorithm, algorithm, level), header=[0])
    x = np.arange(1, 1 + len(df_clusters.columns))
    fig = plt.figure(figsize=(25, 15))
    ax = fig.subplots()
    fraction_bar_plot(x, fraction_sites, ax)
    ax.set_xlabel("cluster", fontsize=35)
    if normalise:
        ax.set_ylabel("fraction of nodes", fontsize=35)
    else:
        ax.set_ylabel("number of nodes", fontsize=35)
    ax.set_title("%s%s distribution across clusters" % ("Shuffled " if shuffled else '', label), fontsize=35)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    n_labels = len(fraction_sites)
    n_col = round(n_labels/20) if n_labels > 20 else 1
    ax.legend(loc='best', bbox_to_anchor=(1, 0.99), fontsize=35, ncol=n_col)
    ax.tick_params(axis='both', labelsize=35)
    plt.show()
    fig.savefig("%s/%s/%s%sclustercomposition_l%d_%s.pdf" % (
        directory, algorithm, "shuffled" if shuffled else '', "fraction_" if normalise else '', int(level), label))


def fraction_bar_plot(x, fraction_sites, ax=None):
    global current_color
    current_color = -1
    if ax is None:
        fig = plt.figure(figsize=(15, 8))
        ax = fig.subplots()
    bottom = np.zeros(len(x))
    ymax = 0
    color_iterator = (color for color in colors_cycle)
    for site, data in fraction_sites.items():
        if np.max(data) == 0:
            continue
        ax.bar(x, data, label=site, bottom=bottom, color=next(color_iterator))
        bottom = bottom + data


def get_Palette(site):
    palette_map = dict({'Brain': 'Blues',
                        'Breast': 'Reds',
                        'Kidney': 'Greens',
                        'Lung': 'Oranges',
                        'Thyroid': 'Greys',
                        'Uterus': 'Purples',
                        'Prostate': 'BuGn',
                        'Ovary': 'BuPu',
                        'Lymph Nodes': 'OrRd',
                        'Soft Tissue': 'PuRd',
                        'Esophagus': 'YlGn',
                        'Stomach': 'YlRd',
                        'Bone Marrow': 'PuBuGn',
                        'Skin': 'YlOrRd',
                        'Adipose Tissue': 'YlOrBr',
                        'Blood': 'RdPu',
                        'Pancreas': 'OrRd',
                        'Testis': 'GnBu'})
    for k in palette_map.keys():
        if k in site:
            return palette_map[k]


current_color = -1


def get_color_cycle():
    global current_color
    current_color += 1
    if current_color >= len(colors_cycle):
        current_color = 0
    return colors_cycle[current_color]


def get_cluster_given_l(l, directory, algorithm='topsbm'):
    df_clusters = pd.read_csv("%s/%s/%s_level_%d_clusters.csv" % (directory, algorithm, algorithm, l), header=[0],
                              index_col=None)
    cluster = {}
    for i, c in enumerate(df_clusters.columns):
        cluster[i] = df_clusters[c].dropna().values
    return cluster


def get_topic_given_l(l, directory, algorithm='topsbm'):
    df_topics = pd.read_csv("%s/%s/%s_level_%d_topics.csv" % (directory, algorithm, algorithm, l), header=[0])
    topic = {}
    for i, c in enumerate(df_topics.columns):
        topic[i] = df_topics[c].dropna().values
    return topic


def get_fraction_sites(cluster, df_files, label='primary_site', normalise=False):
    fraction_sites = {}
    c_fraction_site = {}
    for site in df_files[label].dropna().unique():
        fraction_sites[site] = []
        c_fraction_site[site] = 0

    for i, c in enumerate(cluster):
        for sample in cluster[i]:
            foundsample = get_file(sample, df_files)
            if foundsample is not None:
                c_fraction_site[foundsample[label]] += 1
            else:
                if 'unknown' in fraction_sites.keys():
                    c_fraction_site['unknown'] +=1
                else:
                    c_fraction_site['unknown'] = 1
                    fraction_sites['unknown']=[]
        for site in fraction_sites:
            if normalise:
                norm = float(len(cluster[i]))
            else:
                norm = 1
            if norm > 0:
                fraction_sites[site].append(c_fraction_site[site] / norm)
            else:
                fraction_sites[site].append(np.nan)
            c_fraction_site[site] = 0
    df = pd.DataFrame(data=fraction_sites).dropna(how='all', axis=0)
    ##put first columns that have high values in average
    avgs = df.apply(lambda x: np.average(x.to_numpy()[x.to_numpy().nonzero()[0]]), axis=0)
    df = df.transpose()
    df.insert(0, 'avg', avgs)
    df = df.sort_values(by=['avg'], axis=0, ascending=False).drop('avg', axis=1).transpose()
    df = df.sort_values(by=[tissue for tissue in df.columns], axis=0, ascending=False)
    return df.sort_index(1).to_dict(orient='list')


def get_clustersinfo(cluster, fraction_sites):
    clustersinfo = {
        "maximum": [],
        "homogeneity": [],
        "sizes": [],
        "nclasses": []
    }
    for icluster in cluster:
        maximum = 0
        homo = 0
        size = 0
        nclass = 0
        site_maximum = ''
        cumulative = 0
        for site, data in fraction_sites.items():
            cdata = data[icluster]
            cumulative += cdata
            if cdata > maximum:
                maximum = cdata
                site_maximum = site
            if cdata > 0:
                nclass += 1
                # using fraction_items normalised
                if cdata <= 1:
                    homo -= cdata * np.log(cdata)
            size += cdata
        if cumulative > 0:
            clustersinfo['maximum'].append([float(maximum) / cumulative, site_maximum])
        else:
            clustersinfo['maximum'].append([0, site_maximum])
        clustersinfo['sizes'].append(size)
        clustersinfo['nclasses'].append(nclass)
        clustersinfo['homogeneity'].append(1 - homo)
    return clustersinfo


def plot_maximum(clustersinfo, cluster, label, level, directory, clustersinfo_shuffle=None, algorithm='topsbm'):
    fig = plt.figure(figsize=(15, 6))
    ax = fig.subplots(1, 2)
    bins = 10
    real = np.array(clustersinfo['maximum'])[:, 0].astype(float)
    ax[0].plot(np.sort(real), marker='o', ms=25, ls='')
    ax[1].hist(np.sort(real), histtype='step', bins=bins, lw=4, density=True, range=(0.05, 1.05))
    shuffled = False
    if clustersinfo_shuffle is not None:
        shuffled = np.array(clustersinfo_shuffle['maximum'])[:, 0].astype(float)
        ax[0].plot(np.sort(shuffled), marker='o', ls='', ms=25)
        ax[1].hist(np.sort(shuffled), histtype='step', bins=bins, lw=4, density=True, range=(0.05, 1.05))
        shuffled = True
    ax[0].plot(np.arange(len(cluster)), [0.8 for i in range(len(cluster))], visible=True, ls='--')
    for axi in ax:
        axi.tick_params(axis='both', labelsize=20)
    ax[0].set_xlabel("cluster", fontsize=35)
    ax[0].set_ylabel("maximum fraction\nwith same %s" % label, fontsize=35)
    ax[0].set_ylim((0, 1.1))
    ax[1].set_xlabel("maximum fraction\nwith same %s" % label, fontsize=35)
    ax[1].set_ylabel("pdf", fontsize=35)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    plt.show()
    fig.savefig(
        "%s/%s/%scluster_maximum_l%d_%s.pdf" % (directory, algorithm, "shuffled" if shuffled else '', level, label))


def plot_maximum_size(clustersinfo, label, level, directory, clustersinfo_shuffle=None, algorithm='topsbm'):
    fig = plt.figure(figsize=(15, 6))
    x = np.array(clustersinfo['sizes']).astype(int)
    y = np.array(clustersinfo['maximum'])[:, 0].astype(float)
    plt.scatter(x, y, lw=10, label='clusters')
    plt.xlim(0, np.max(x) + np.max(x) / 10)
    plt.plot(np.linspace(0.5, x.max()), 1. / np.linspace(0.5, x.max()), label='uniform')
    shuffled = False
    if clustersinfo_shuffle is not None:
        shuffled = True
        x_shuffle = np.array(clustersinfo_shuffle['sizes']).astype(int)
        y_shuffle = np.array(clustersinfo_shuffle['maximum'])[:, 0].astype(float)
        plt.scatter(x_shuffle, y_shuffle, lw=10, label='clusters shuffled')
        plt.xlim(0, np.max(x_shuffle) + np.max(x_shuffle) / 10)
    plt.xlabel("cluster size", fontsize=35)
    plt.ylabel("maximum fraction\nwith same %s" % label, fontsize=35)
    plt.ylim((0, 1.1))
    plt.legend(loc='best', fontsize=35)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    plt.show()
    fig.savefig(
        "%s/%s/%sclusterhomosize_l%d_%s.pdf" % (directory, algorithm, "shuffled" if shuffled else '', level, label))


def plot_maximum_label(clustersinfo, label, level, directory, clustersinfo_shuffle=None, algorithm='topsbm'):
    fig = plt.figure(figsize=(10, 6))
    x = np.array(clustersinfo['nclasses']).astype(int)
    y = np.array(clustersinfo['maximum'])[:, 0].astype(float)
    shuffled = False
    plt.scatter(x, y, lw=10, alpha=0.9, label='clusters')
    plt.plot(np.arange(1, np.max(x) + 2), 1. / np.arange(1, np.max(x) + 2), ls='--', c='cyan', label='uniform')
    plt.xlim(0.95, np.max(x) + 0.5)
    if clustersinfo_shuffle is not None:
        x_shuffle = np.array(clustersinfo_shuffle['nclasses']).astype(int)
        y_shuffle = np.array(clustersinfo_shuffle['maximum'])[:, 0].astype(float)
        plt.scatter(x_shuffle, y_shuffle, lw=10, alpha=0.9, label='clusters shuffled')
        plt.plot(np.arange(1, np.max(x_shuffle) + 2), 1. / np.arange(1, np.max(x_shuffle) + 2), ls='--', c='cyan',
                 label='')
        shuffled = True
        plt.xlim(0.95, np.max(x_shuffle) + 0.5)
    plt.xlabel("number of labels", fontsize=35)
    plt.ylabel("maximum fraction\nwith same %s" % label, fontsize=35)
    plt.ylim((0, 1.1))
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.legend(loc='lower right', fontsize=35)
    plt.show()
    fig.savefig(
        "%s/%s/%scluster_homon_l%d_%s.pdf" % (directory, algorithm, "shuffled" if shuffled else '', level, label))


def plot_labels_size(clustersinfo, label, level, directory, clustersinfo_shuffle=None, algorithm='topsbm'):
    fig = plt.figure(figsize=(10, 6))
    x = np.array(clustersinfo['sizes']).astype(float)
    y = np.array(clustersinfo['nclasses']).astype(int)
    plt.xlim(x.min() - 10, x.max() + 5)
    plt.ylim(y.min() - 2, y.max() + 5)
    shuffled = False
    plt.scatter(x, y, lw=10, alpha=0.9, label='clusters')
    if clustersinfo_shuffle is not None:
        x_shuffle = np.array(clustersinfo_shuffle['sizes']).astype(float)
        y_shuffle = np.array(clustersinfo_shuffle['nclasses']).astype(int)
        plt.scatter(x_shuffle, y_shuffle, lw=10, alpha=0.9, label='clusters shuffled')
        plt.xlim(x.min() - 10, x_shuffle.max() + 5)
        plt.ylim(y.min() - 2, y_shuffle.max() + 8)
        shuffled = True
    plt.xlabel("cluster size", fontsize=35)
    plt.ylabel("number of labels", fontsize=35)
    plt.legend(loc='upper right', fontsize=35)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.show()
    fig.savefig(
        "%s/%s/%scluster_shuffle_label_size_l%d_%s.pdf" % (
            directory, algorithm, "shuffled" if shuffled else '', level, label))


def make_heatmap(fraction_sites, directory, label, level, shuffled=False, normalise=False, algorithm='topsbm'):
    sns.set(font_scale=2)
    found_classes = []
    for site, data in fraction_sites.items():
        if np.max(data) == 0:
            continue
        found_classes.append(site)
    for arr in fraction_sites.values():
        x = len(arr)
        break
    x = np.arange(1, 1 + x)
    fig = plt.figure(figsize=(30, 10))
    fig.subplots(1)
    sns.heatmap(pd.DataFrame(data=fraction_sites).loc[:, found_classes].transpose(), vmin=0, cmap="RdYlBu_r",
                xticklabels=x)
    fig.savefig("%s/%s/%sheatmap_cluster%s_l%d_%s.pdf" % (
        directory, algorithm, "shuffled" if shuffled else '', "fraction_" if normalise else '', int(level), label))


def get_file(sample, df_file):
    for fullsample in df_file.index.values:
        if sample in fullsample:
            return df_file.loc[fullsample, :]
    return None


def define_labels(cluster, df_files, label='primary_site', verbose=False):
    true_labels = []
    predicted_labels = []
    for c in cluster:
        if verbose:
            print(c)
        for sample in cluster[c]:
            try:
                true_labels.append(get_file(sample, df_files)[label])
                predicted_labels.append(c)
            except:
                true_labels.append('')
                predicted_labels.append('')
                print(*sys.exc_info())
                print("error searching %s in %s" % (label, sample))
    _, true_labels = np.unique(true_labels, return_inverse=True)
    return true_labels, predicted_labels


def add_score_lines(ax, scores, V="V", labels=None, h=False, c=False, alpha=0.8, **kwargs):
    '''
    add to ax lines in scores
    add homogeneity and completness if required by h and c
    '''
    colors = {
        'primary_site': 'blue',
        'hsbm': 'blue',
        'secondary_site': 'red',
        'status': 'red',
        'hSBM': 'blue',
        'mixed': 'green',
        'hierhsbm': 'purple',
        'hsbm->hierachical': 'purple',
        'disease_type': 'red',
        'shuffle': 'orange',
        'tm': 'darkcyan',
        'cc': 'darkred',
        'disease_tissue': 'purple',
        'hierarchical': 'darkgreen',
        'lda': 'violet',
        'RPPA Clusters': 'red',
        'wgcna': 'purple',
        "Subtype_Selected": "red",
        "BRCA_Subtype_PAM50": "blue"
    }

    for label in labels:
        if label not in scores.keys():
            print("No score for %s"%label)
            continue
        if label not in colors.keys():
            colors[label]=next(color_iterator)
        xl = scores[label]['xl']
        if h:
            ax.plot(xl, scores[label]['h'], ls='-.', c=colors[label], marker='x', lw=150, ms=45, alpha=alpha,
                    label='homogeneity - %s' % label)
        if c:
            ax.plot(xl, scores[label]['c'], ls=':', c=colors[label], marker='<', lw=10, ms=45, alpha=alpha,
                    label='completness - %s' % label)
        if len(scores[label][V]) == len(xl):
            ax.plot(xl, scores[label][V], label='%s' % label, ls='-', c=colors[label], marker='o', lw=20, ms=45,
                    **kwargs)
        else:
            raise(ValueError("xl has got wrong lenght"))
    customize_metric_plot(ax, xl)


def customize_metric_plot(ax, xl):
    ax.tick_params(labelsize=35, width=8, length=20)
    ax.tick_params(which="minor", labelsize=35, width=5, length=15)
    ax.set_xlabel("Number of clusters", fontsize=40)
    ax.set_ylabel("NMI score", fontsize=40)
    ax.set_ylim((0, 1.1))
    ax.set_xlim(1, np.max(xl)*1.1)
    ax.set_xscale('log')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='best', bbox_to_anchor=(1, 0.85), fontsize=35, ncol=1)

def plot_topic_size(directory, l, algorithm='topsbm'):
    df_topics = pd.read_csv("%s/%s/%s_level_%d_topics.csv" % (directory, algorithm, algorithm, l))
    sizes = []
    for t in df_topics.columns:
        sizes.append(len(df_topics.loc[:, t].dropna()))
    bins = np.linspace(0.5, np.max(sizes) + 0.5, int((np.max(sizes) + 1) / (np.max(sizes) / 5)))
    bin_counts, bin_edges, _ = plt.hist(sizes, histtype='step', lw=2, bins=bins)
    fig = plt.figure()
    ax = fig.subplots()
    ax.set_title("[%d topics, level: %d]" % (len(df_topics.columns), l))
    x = (bin_edges[:-1] + bin_edges[1:]) / 2
    ax.plot(x[np.nonzero(bin_counts)], bin_counts[np.nonzero(bin_counts)])
    ax.plot(x, 1e4 / np.power(x, 5))
    ax.set_xlabel("topic size\n(number of genes)", fontsize=35)
    ax.set_ylabel("number of topic", fontsize=35)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.show()
    fig.savefig("%s/%s/topic_size_level%d.png" % (directory, algorithm, l))


def get_candles(directory, level, df_mv, ax, algorithm='topsbm'):
    df_topics = pd.read_csv("%s/%s/%s_level_%d_topics.csv" % (directory, algorithm, algorithm, level))
    candles = {
        'open': [],
        'high': [],
        'low': [],
        'close': [],
        'size': []
    }
    for topic in df_topics.columns:
        subarr = df_mv.loc[df_topics[topic].dropna(), :]['occurrence'].values
        avg = np.average(subarr)
        std = np.std(subarr)
        q = np.quantile(subarr, [0.25, 0.75])
        candles['high'].append(np.min([1, avg + std]))
        candles['open'].append(np.min([q[1], 1]))
        candles['close'].append(np.max([q[0], 0]))
        candles['low'].append(np.max([0, avg - std]))
        candles['size'].append(len(subarr))
    ax.set_title("[level: %d]" % level)
    ax.set_ylabel('$O_i$', fontsize=35)
    ax.set_xlim(-1, len(df_topics.columns))
    ax.set_xticks([i + 1 for i in range(-1, len(df_topics.columns))])
    ax.set_xticklabels(
        ["Topic %d" % (i + 2) if ((i + 2) % 5 == 0 or i == -1) else '' for i in range(-1, len(df_topics.columns))],
        rotation=60)
    return candles


def get_tissue_style(tissue):
    marker = 'o'
    c = 'k'
    ls = '--'
    if 'gtex' in tissue:
        marker = 'o'
        ls = '-'
    elif 'tcga' in tissue:
        marker = 'x'
        ls = '--'
    else:
        marker = '.'
        ls = '-.'
    if 'reast' in tissue:
        c = 'darkcyan'
    elif 'olon' in tissue:
        c = 'b'
    elif 'hyroid' in tissue:
        c = 'y'
    elif 'terus' in tissue:
        c = 'pink'
    elif 'ladder' in tissue:
        c = 'gray'
    elif 'sophagus' in tissue:
        c = 'brown'
    elif 'ung' in tissue:
        c = 'magenta'
    elif 'tomach' in tissue:
        c = 'lime'
    elif 'kin' in tissue:
        c = 'wheat'
    elif 'ancreas' in tissue:
        c = 'forestgreen'
    elif 'Adrenal Gland' in tissue:
        c = 'aqua'
    elif 'Adipose Tissue' in tissue:
        c = 'brown'
    elif 'erve' in tissue:
        c = 'royalblue'
    elif 'lood' in tissue:
        c = 'red'
    elif 'idney' in tissue:
        c = 'mediumslateblue'
    elif 'eart' in tissue:
        c = 'darkred'
    elif 'rain' in tissue:
        c = 'darkgray'
    elif 'estis' in tissue:
        c = 'darkkhaki'
    elif 'LumA' in tissue:
        c = 'pink'
    elif 'LumB' in tissue:
        c = 'purple'
    elif 'Normal' in tissue:
        c = 'blue'
    elif 'Basal' in tissue:
        c = 'darkred'
    elif 'Her2' in tissue:
        c = "green"
    else:
        c = 'k'
    return (marker, c, ls)


def topic_distr_sample(doc, df, ax=None):
    if ax == None:
        fig = plt.figure()
        ax = fig.subplots()
    ax.set_title("Topic distribution: %s" % doc)
    labels = [l if df[df['doc'] == doc].loc[:, l].values[0] >= 0.05 else '' for l in df.columns[2:]]
    patches, texts, autotexts = ax.pie(df[df['doc'] == doc].values[0][2:], labels=labels,
                                       autopct=lambda p: '%.1f%s' % (p, '%') if p >= 5 else '',
                                       textprops={'fontsize': 20, 'color': 'white', 'wrap': True})
    for t in texts:
        t.set_fontsize(18)
        t.set_wrap(True)
        t.set_color('black')
    plt.show()


def topic_distr_isample(idoc, df, ax=None):
    topic_distr_sample(df[df['i_doc'] == idoc]['doc'].values[0], ax)


def add_tumor_location(df_files):
    df_files.insert(2, 'disease_tissue', '')
    for sample in df_files.index.values:
        row = df_files.loc[sample, :]
        df_files.at[sample, 'disease_tissue'] = '%s[%s]' % (row['primary_site'], row['disease_type'])


def get_scores(directory, labels, df_files=None, algorithm='topsbm', verbose=False):
    if df_files is None:
        df_files = pd.read_csv("%s/files.dat" % directory, index_col=[0], header=[0]).dropna(how='all', axis=0)
    if df_files.columns.isin(['disease_type']).any():
        add_tumor_location(df_files)
    scores = {}
    for label in labels:
        xl = []
        scores[label] = {
            'h': [],
            'c': [],
            'V': [],
            'xl':[]
        }
        l = get_max_available_L(directory, algorithm)
        for l in np.arange(l + 1):
            try:
                true_labels, predicted_labels = define_labels(get_cluster_given_l(l, directory, algorithm), df_files, label=label)
                scores[label]['h'].append(metrics.cluster.homogeneity_score(true_labels, predicted_labels))
                scores[label]['c'].append(metrics.cluster.completeness_score(true_labels, predicted_labels))
                scores[label]['V'].append(metrics.cluster.v_measure_score(true_labels, predicted_labels))
                xl.append(len(np.unique(predicted_labels)))
                if verbose:
                    print(l)
            except:
                print(*sys.exc_info())
                print("Skipping level ", l)

        # add the first point where all sample are in the same cluster by definition
        if 1 not in xl:
            if xl[0] < xl[-1]:
                idx = 0
            else:
                idx = len(xl)
            true_labels, _ = define_labels(get_cluster_given_l(l, directory, algorithm), df_files, label=label)
            predicted_labels = np.ones_like(true_labels)
            scores[label]['h'].insert(idx,metrics.cluster.homogeneity_score(true_labels, predicted_labels))
            scores[label]['c'].insert(idx,metrics.cluster.completeness_score(true_labels, predicted_labels))
            scores[label]['V'].insert(idx,metrics.cluster.v_measure_score(true_labels, predicted_labels))
            xl.insert(idx,len(np.unique(predicted_labels)))

        scores[label]['xl'] = xl
    if len(labels) >= 2:
        h = np.array(scores[labels[0]]['h'])
        c = np.array(scores[labels[1]]['c'])
        scores['mixed'] = {
            'h': h,
            'c': c,
            'V': 2 * h * c / (h + c)
        }
    return scores

def shuffle_files(df_files, label, random_state=42):
    df_files_shuffled = df_files.copy()
    if label not in df_files.columns:
        raise(AttributeError(f"{label} non available in:{df_files.columns}"))
    df_files_shuffled[label] = shuffle(df_files_shuffled[label].values)
    return df_files_shuffled


def get_scores_shuffled(directory, df_files, algorithm='topsbm', label='primary_site', verbose=False):
    scores = {
        'h': [],
        'c': [],
        'V': [],
        'xl':[]
    }
    xl = []
    l = get_max_available_L(directory, algorithm)
    df_files_shuffled = shuffle_files(df_files.copy(), label, random_state=42)
    try:
        for l in np.arange(0, l + 1):
            try:
                if verbose:
                    print(l)
                clusters = get_cluster_given_l(l, directory, algorithm=algorithm)
            except:
                print("Skipping shuffled level ", l)
                continue
            _, predicted_labels = define_labels(clusters, df_files, label=label)
            true_labels, _ = define_labels(clusters,
                                           df_files_shuffled,
                                           label=label)
            scores['h'].append(metrics.cluster.homogeneity_score(true_labels, predicted_labels))
            scores['c'].append(metrics.cluster.completeness_score(true_labels, predicted_labels))
            scores['V'].append(metrics.cluster.v_measure_score(true_labels, predicted_labels))
            xl.append(len(np.unique(predicted_labels)))
    except:
        print(*sys.exc_info())
        print("shuffled files not found")

    # add the first point where all sample are in the same cluster by definition
    if xl[0] < xl[-1]:
        idx = 0
    else:
        idx = len(xl)
    true_labels, _ = define_labels(get_cluster_given_l(l, directory, algorithm), df_files, label=label)
    predicted_labels = np.ones_like(true_labels)
    scores['h'].insert(idx,metrics.cluster.homogeneity_score(true_labels, predicted_labels))
    scores['c'].insert(idx,metrics.cluster.completeness_score(true_labels, predicted_labels))
    scores['V'].insert(idx,metrics.cluster.v_measure_score(true_labels, predicted_labels))
    xl.insert(idx,len(np.unique(predicted_labels)))

    scores['xl'] = xl
    return scores


def getclustersizesarray(directory, l=3, algorithm='topsbm'):
    try:
        xl = [len(get_cluster_given_l(li, directory, algorithm=algorithm)) for li in np.linspace(0, l, l + 1)]
    except:
        try:
            xl = [len(get_cluster_given_l(li, directory, algorithm=algorithm)) for li in np.linspace(1, l, l)]
        except:
            xl = []
            for li in np.linspace(1, l, l):
                try:
                    xl.append(len(get_cluster_given_l(li, directory, algorithm=algorithm)))
                except:
                    pass
    return xl


def gettopicsizesarray(directory, l=3, algorithm='topsbm'):
    xl = []
    try:
        xl = [len(get_topic_given_l(li, directory, algorithm=algorithm)) for li in np.linspace(0, l, l + 1)]
    except:
        try:
            xl = [len(get_topic_given_l(li, directory, algorithm=algorithm)) for li in np.linspace(1, l, l)]
        except:
            xl = []
            for li in np.linspace(1, l, l):
                try:
                    xl.append(len(get_topic_given_l(li, directory, algorithm=algorithm)))
                except:
                    pass
    return xl


def plot_sizes(level, directory, algorithm, ax=None):
    cluster = get_cluster_given_l(level, directory, algorithm=algorithm)
    if ax is None:
        fig = plt.figure(figsize=(10, 6))
        ax = fig.subplots()
    sizes = []
    for c in cluster.items():
        sizes.append(len(c[1]))
    ax.set_xlabel("size", fontsize=24)
    ax.set_ylabel("number of clusters", fontsize=24)
    ax.set_title("Cluster sizes at level %d" % level)
    ax.hist(sizes, histtype='step', lw=4)
    plt.savefig("%s/%s/sizes_distr_level%d.pdf" % (directory, algorithm, level))
    plt.show()


def clusteranalysis(directory, labels, algorithm='topsbm') -> None:
    """
    Perform analyses of an algorithm output

    :param directory: where to search the data
    :param labels: ground truth label to search. This should be in a file called directory/files.dat
    :param algorithm: name of the folder in which data are stored
    """
    l_max = get_max_available_L(directory, algorithm)
    df_clusters = pd.read_csv("%s/%s/%s_level_%d_clusters.csv" % (directory, algorithm, algorithm, l_max), header=[0])
    if df_clusters is None:
        print("files not found")
    df_files = pd.read_csv("%s/files.dat"%directory, index_col=[0], header=[0]).dropna(axis=1, how='all').dropna(axis=0, how='all')
    samples = pd.read_csv("%s/%s/%s_level_0_clusters.csv"%(directory,algorithm,algorithm), header=[0]).astype(str).values.ravel()
    samples=samples[samples!="nan"]
    df_files = df_files.reindex(index=samples).dropna(how="all", axis=0).fillna("unknown")
    df_files_shuffled = df_files.copy()
    df_files_shuffled.apply(lambda x: np.random.shuffle(x), 0)
    for normalise in [True, False]:
        for label in labels:
            for level in np.arange(l_max+1)[::-1]:
                print(normalise, label, level)
                try:
                    cluster = get_cluster_given_l(level, directory,algorithm=algorithm)
                    fraction_sites = get_fraction_sites(cluster,df_files=df_files,label=label, normalise=normalise)

                    clustersinfo = get_clustersinfo(cluster,fraction_sites)
                    plot_cluster_composition(fraction_sites,directory,level,label=label, normalise=normalise,algorithm=algorithm)
                    make_heatmap(fraction_sites, directory, label, level, normalise=normalise,algorithm=algorithm)

                    clustersinfo = get_clustersinfo(cluster,fraction_sites)
                    if not normalise:
                        plot_maximum(clustersinfo,cluster,label,level, directory,algorithm=algorithm)
                        plot_maximum_size(clustersinfo,label,level, directory,algorithm=algorithm)
                        plot_maximum_label(clustersinfo,label,level, directory,algorithm=algorithm)
                        plot_sizes(level,directory, algorithm=algorithm)
                except:
                    print(*sys.exc_info())
                continue
                shuffle_files(df_files,label).to_csv("%s/files_shuffles.dat"%directory, index=True)
                fraction_sites_shuffle = get_fraction_sites(cluster, df_files=pd.read_csv("%s/files_shuffles.dat"%directory, index_col=[0]),label=label, normalise=normalise)
                clustersinfo_shuffle = get_clustersinfo(cluster, fraction_sites_shuffle)
                plot_cluster_composition(fraction_sites_shuffle,directory,level, label=label, shuffled=True, normalise=normalise, algorithm=algorithm)
                if not normalise:
                    plot_maximum(clustersinfo,cluster,label,level,directory,clustersinfo_shuffle,algorithm=algorithm)
                    plot_maximum_size(clustersinfo,label,level, directory,clustersinfo_shuffle,algorithm=algorithm)
                    plot_maximum_label(clustersinfo,label,level, directory,clustersinfo_shuffle,algorithm=algorithm)
                    plot_labels_size(clustersinfo,label,level, directory,clustersinfo_shuffle,algorithm=algorithm)
    ##define scores
    scores = get_scores(directory, labels, algorithm=algorithm)
    try:
        xl = getclustersizesarray(directory, l_max)
        with open("%s/clustersizes.txt" % directory, 'w') as f:
            for x in xl:
                f.write("%d\n" % x)
    except:
        print("cannot save clustersizes.txt")

    try:
        xl = gettopicsizesarray(directory, l_max)
        with open("%s/topicsizes.txt" % directory, 'w') as f:
            for x in xl:
                f.write("%d\n" % x)
    except:
        print("cannot save topicsizes.txt")

    # save files for R analisys
    for l_max in np.arange(l_max + 1):
        pd.DataFrame(data=define_labels(get_cluster_given_l(l_max, directory, algorithm=algorithm), df_files, label=labels[0])[1],
                     columns=['l%d' % l_max]).to_csv("%s/%s/%s_level_%d_labels.csv" % (directory, algorithm, algorithm, l_max),
                                                     header=True,
                                                     index=False)


def get_max_available_L(directory, algorithm='topsbm'):
    """
    Get maximum layer available for algorithm
    """
    return np.array([el.split("_")[2] for el in os.listdir("%s/%s" % (directory, algorithm)) if "level_" in el],
                    dtype=int).max()

def out_to_file(out, index, name='new_method', l=0):
    print("saving clusters")
    df_clusters = pd.DataFrame(index=np.arange(len(index)))
    for c in np.arange(out.max()+1)[::-1]:
        c_objects = index[np.argwhere(out == c)].values.T[0]
        df_clusters.insert(0, "Cluster %d" % (c + 1),
                           np.concatenate((c_objects, [np.nan for _ in np.arange(len(index) - len(c_objects))])))
    df_clusters.dropna(axis=0, how='all', inplace=True)
    df_clusters.to_csv("%s_level_%d_clusters.csv"%(name, l), index=False, header=True)


#normalise to hsbm
def normalise_score(scores : dict, base_algorithm="hsbm", operation=lambda x,y: x/y, epsilon = 1e-6)->None:
    "save scaled data to scores[norm_V]"
    for algorithm in scores.keys(): #the first point is always constructed and np.interp wants sorted data so[:-1:-1]
        baseline = np.interp(scores[algorithm]["xl"],
    					 scores[base_algorithm]["xl"][:-1][::-1],
    					 scores[base_algorithm]["V"][:-1][::-1])
        scores[algorithm]["norm_V"]=operation(np.array(scores[algorithm]["V"])+epsilon,baseline+epsilon)
