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

import sys

import gseapy as gs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from topicpy.tableanalyser import get_symbol


def get_ontology_df(topic, cutoff=0.05, threshhold=5e-1, gene_sets = ['GO_Molecular_Function_2018',
             'GO_Biological_Process_2018',
             'GO_Cellular_Component_2018',
             'Human_Phenotype_Ontology',
             'GTEx_Tissue_Sample_Gene_Expression_Profiles_up',
             'GTEx_Tissue_Sample_Gene_Expression_Profiles_down',
             'Tissue_Protein_Expression_from_Human_Proteome_Map',
             'KEGG_2019_Human',
             'NCI-60_Cancer_Cell_Lines'
            ], background=None) -> pd.DataFrame:
    """

    :param topic: list of genes
    :param background: enrichment test background
    :param cutoff: Enrichments cutoff
    :param threshhold: threshold on Adjusted P-value
    :return: DataFrame with terms and P-vals
    """
    sets = ','.join(gene_sets)
    if background is None:
        background='hsapiens_gene_ensembl'
    topic = [g for g in topic if str(g)!='nan']
    gene_ontology = gs.enrichr(list(topic), gene_sets=sets, cutoff=cutoff, background=background).results
    return gene_ontology[gene_ontology['Adjusted P-value'] < threshhold][['Term', 'Adjusted P-value', 'Gene_set']]


def ensg_to_symbol(series):
    """
    convert gene ENSG to Symbol
    :param series: series of ENSGs
    :return: symbols
    """
    symbols = []
    for g in series:
        try:
            symbols.append(get_symbol(g))
        except:
            print(*sys.exc_info())
    symbols = np.array(symbols)
    return symbols[[s!='nan' for s in symbols]]


def save_plot_Pvalues(df_topics, l, directory, algorithm):
    """

    :param df_topics: Topics DataFrame
    :param l: level of the analysis
    :param directory:
    """
    topic_pvalues = []
    topic_gos = []
    for itopic, topic in enumerate(df_topics.columns):
        try:
            enriched_topic = pd.read_csv("%s/%s/gsea_level_%d_topic_%d.csv" % (directory, algorithm, l, itopic+1))
            if len(enriched_topic.index) > 0:
                p_val = np.sort(enriched_topic['Adjusted P-value'])[0]
                topic_pvalues.append(-np.log10(p_val))
                for goc in enriched_topic['Gene_set'][:10].unique():
                    topic_gos.append(goc)
        except:
            print(*sys.exc_info())

    fig = plt.figure()
    #x = np.arange(1, 1 + len(topic_pvalues))
    c, _, _ = plt.hist(topic_pvalues, histtype='step', lw=2)
    plt.plot([-np.log10(0.05) for _ in np.linspace(1, 10, num=10)], np.arange(0, np.max(c) + 5, (np.max(c) + 5) / 10),
             ls='--', lw=5, label="$\\alpha=0.05$")
    plt.xlabel('-log(P-value)', fontsize=20)
    plt.ylabel("number of topics", fontsize=20)
    # plt.ylim(0,0.055)
    # plt.yscale('log')
    plt.legend(fontsize=18)
    fig.savefig("%s/pvaluescrosstopic(%d).pdf" % (directory, l))

    fig = plt.figure(figsize=(20, 10))
    gos, goscounts = np.unique(topic_gos, return_counts=True)
    from textwrap import wrap
    plt.barh(["\n".join(wrap(str(l).replace('_', ' '), 20)) for l in gos], goscounts)
    plt.yticks(fontsize=15)
    plt.show()
    fig.savefig("%s/pvaluecategories(%d).pdf" % (directory, l))




def topic_analysis(directory, l, algorithm="topsbm", background=None, verbose=True, save_Pvalues=True):
    """

    :param directory: where to find files
    :param l: level of the analisys
    :param algorithm: name of folder and files containing topics table (e.g. path/to/topsbm/topsbm_level_3_topics.csv)
    :param background: List of background genes
    :param verbose: verbosity
    :param save_Pvalues: save data for P-values plot
    """
    df_topics = pd.read_csv("%s/%s/%s_level_%d_topics.csv" % (directory, algorithm, algorithm, l))
    enriched_topic = None
    for itopic, topic in enumerate(df_topics.columns):
        try:
            enriched_topic = pd.read_csv("%s/%s/gsea_level_%d_topic_%d.csv" % (directory, algorithm, l, itopic+1), index_col=[0])
            print(topic)
        except:
            try:
                symbols = ensg_to_symbol(df_topics.loc[:, topic].dropna().values)
                print(topic)
                enriched_topic = get_ontology_df(symbols, background=background).sort_values(by=['Adjusted P-value'],
                                                                                             ascending=True)
                enriched_topic = enriched_topic.loc[enriched_topic.index.values[:20], :]
                enriched_topic.to_csv("%s/%s/gsea_level_%d_topic_%d.csv" % (directory, algorithm, l, itopic+1))
            except:
                print(*sys.exc_info())
        if verbose:
            print(enriched_topic)
    if save_Pvalues:
        save_plot_Pvalues(df_topics, l, directory, algorithm)
