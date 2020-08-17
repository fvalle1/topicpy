import pandas as pd
import numpy as np
from scipy.stats import hypergeom 


def parameters_for_hypergeometric(list_1, list_2):
    """
    Returns
    x num of successes
    M population size
    k successes in population
    N sample size
    """
    population_size = len(list_1[list_1.index.isin(list_2.index)])
    pop_successes = {module:len(list_2[list_2==module]) for module in list_2.unique()}
    sample_sizes = {topic:len(list_1[list_1==topic]) for topic in list_1.unique()}
    num_successes = pd.DataFrame(index=list_1.unique(), columns=list_2.unique()).fillna(0)
    for g in list_2.index:
        if g in list_1.index:
            num_successes.at[list_1[g],list_2[g]]+=1
    print(num_successes.shape)
    return num_successes, population_size, pop_successes, sample_sizes, (list_1, list_2)

def build_map(num_successes, population_size, pop_successes, sample_sizes, lists, last_name=None):
    list_1, list_2 = lists
    df_cmap=pd.DataFrame(index=["Topic %d"%(d+1) for d in range(len(list_1.unique()))],
                     columns=["Topic %d"%(d+1) for d in range(len(list_2.unique()))]).fillna(0.5)
    if last_name is not None:
        df_cmap.columns=num_successes.columns
    for module, module_successes in zip(df_cmap.columns, num_successes.columns):
        for topic, topic_successes in zip(df_cmap.index, num_successes.index):
            x = num_successes.at[topic_successes,module_successes].astype(int) # number of successes
            M = population_size # pop size
            k = pop_successes[module_successes] # successes in pop
            N = sample_sizes[topic_successes] # sample size
            pval = hypergeom.sf(x-1, M, k, N)
            df_cmap.at[topic,module]=-np.log10(float(pval))
    return df_cmap


def plot_map(df_cmap, first_name="topsbm", last_name="lda", *args, **kwargs):
    import seaborn as sns
    import matplotlib.pyplot as plt

    #df_cmap = df_cmap.sort_values(by=[c for c in df_cmap.columns], axis=0, ascending=True)
    #create a color palette with the same number of colors as unique values in the Source column
    network_pal = sns.color_palette('husl',n_colors=len(df_cmap.columns))

    #Create a dictionary where the key is the category and the values are the
    #colors from the palette we just created
    network_lut = dict(zip(df_cmap.columns, network_pal))
    network_col = df_cmap.columns.map(network_lut)
    #Create a dictionary where the key is the category and the values are the
    #colors from the palette we just created
    network_lut = dict(zip(df_cmap.columns, network_pal))
    network_col = df_cmap.columns.map(network_lut)

    cm = sns.clustermap(df_cmap, 
                        row_cluster=False, 
                        col_cluster=False, 
                        metric='euclidean', 
                        vmin=0,
                        vmax = 30,
                        cmap='Blues_r', 
                        col_colors=network_col,
                        mask=False,
                        cbar_pos=(1.05,0.05,0.05,0.7),
                        *args, 
                        **kwargs)

    ax = cm.ax_heatmap
    ax.tick_params(labelsize=15)

    ax.set_ylabel(first_name, fontsize=35)
    ax.set_xlabel(last_name, fontsize=35)
    ax.set_xticklabels(["Topic %d"%(t+1) for t,_ in enumerate(df_cmap.columns)], rotation=75)
    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")
    ax.set_yticklabels(["Topic %d"%(t+1) for t,_ in enumerate(df_cmap.index)], rotation=0)

    cax = cm.ax_cbar
    cax.tick_params(labelsize=35)
    cax.set_title("-Log(P-value)", fontsize=30)

    #plt.tight_layout()

    cm.fig.suptitle('Algorithm comparison', fontsize=40)
    cm.savefig(f"topics_logp_{first_name}_{last_name}.pdf")
    plt.show()
    return cm