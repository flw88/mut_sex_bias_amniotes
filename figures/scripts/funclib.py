#!/usr/bin/env python

import pandas as pd
import numpy as np
import itertools
import sys
import matplotlib.pyplot as plt
import statsmodels.api as sm
import math
import statsmodels.stats.proportion as ssp
import matplotlib.colors as colors_mplt
from adjustText import adjust_text
import statsmodels.stats.proportion

def df2dict(df, key, value):
    '''Returns dictionary out of two columns in a dataframe'''
    return df.set_index(key)[value].to_dict()

def min_div_to_chrom_assembly(df, target_sp, clevel_sps, tree):
    '''Returns the minimum distance between sp and a chromosome level assembly'''
    if target_sp in clevel_sps:
        return 0
    distances = [tree.get_distance(target_sp, sp) for sp in clevel_sps]
    return min(distances)

def get_worst_species(df, decision_tree):
    '''Returns species ranking higher in decision tree'''
    for param, criteria in decision_tree.items():
        if criteria!="Max" and criteria!="Min":
            subd = df[df[param]!=criteria]
            if len(subd)==1:
                return subd.Species.values[0]
        else:
            if len(df[param].unique())==1:
                continue
            sorted_d = df.sort_values(by=param, ascending=False if criteria!="Max" else True)
            sorted_d.index = range(len(sorted_d))
            return sorted_d.Species.values[0]

def conservative_prune_dataset_by_pi(df, times_pi, decision_tree, tree):
    '''Thins dataset, ensuring species are at list times_pi*pi diverged and >2%
    keeping species according to decision tree'''
    
    exclude_species = []
    sp_universe = [tree.get_leaf_names()]
    for sp1,sp2 in itertools.combinations(df.Species, 2):
        distance = tree.get_distance(sp1,sp2)
        pi1 = float(df[df.Species==sp1]["Pi_het"])
        pi2 = float(df[df.Species==sp2]["Pi_het"])
        max_pi = np.nanmax([pi1,pi2])
        min_distance = max_pi*times_pi
        if min_distance>=distance or distance<=0.02:
            worst_sp = get_worst_species(df[df.Species.isin([sp1,sp2])], decision_tree)
            exclude_species.append(worst_sp)
        
    pruned_df = df[~df.Species.isin(exclude_species)]
    return pruned_df
    
def prune_dataset_by_pi(df, times_pi, decision_tree, tree):
    '''Thins dataset, ensuring species are at list times_pi*pi diverged
    keeping species according to decision tree'''
    
    # All comparisons, record if pair is valid
    total_species = df.Species
    matrix = []
    for sp1,sp2 in itertools.combinations(total_species, 2):
        distance = tree.get_distance(sp1,sp2)
        pi1 = float(df[df.Species==sp1]["Pi_het"])
        pi2 = float(df[df.Species==sp2]["Pi_het"])
        max_pi = np.nanmax([pi1,pi2]) if not all(math.isnan(p) for p in [pi1,pi2]) else 10
        min_distance = max_pi*times_pi
        matrix.append([sp1,sp2,0 if min_distance>distance else 1, distance])
        
    c_df = pd.DataFrame(matrix)
    c_df.columns = ["sp1","sp2","acceptance","divergence"]
    return c_df
    # Find valid combination with highest numerber of species
    for n in range(len(total_species)+1)[::-1]:
        sys.stderr.write("Trying combination of {} species...\n".format(n))
        for sp_set in itertools.combinations(total_species, n):
            set_df = c_df[(c_df.sp1.isin(sp_set)) & 
                          (c_df.sp2.isin(sp_set))]
            validity = set_df.acceptance.sum()==len(set_df)
            if validity:
                sys.stderr.write("Found a valid set of {} species.\n".format(n))
                return sp_set

    sys.stderr.write("There's no safe dataset...!")

def replace_values(array,minv,maxv):
    '''Truncate range of values given a min and max value'''
    new_array = []
    for a in array:
        if a<=minv:
            new_array.append(minv)
        elif a>=maxv:
            new_array.append(maxv)
        else:
            new_array.append(a)
    return new_array

def add_columns_to_df(og_df, add_df, join_col, cols):
    '''Adds columns (col) from a df (add_df) to another (og_df) 
    based on a shared variable (join_col)'''

    joincol_to_col = add_df.set_index(join_col)[cols].to_dict()
    
    for c in cols:
        og_df[c] = [joincol_to_col[c][sp] if sp in joincol_to_col[c] else np.nan for sp in og_df[join_col]]

    return og_df

def miyata(xar):
    """Alpha estimate from X/A subsitution rate ratio"""
    return (4-(3*xar))/((3*xar)-2)

def reverse_miyata(alpha):
    """Given alpha, return X/A ratio"""
    return (2*alpha+4)/float(3*alpha+3)

def alpha_from_pedigrees(dnm_file):
    """Reads dnms from file and computes alpha and binomial CIs"""
    dnm_df = pd.read_csv(dnm_file,sep="\t")
    dnm_df["Species_simple"] = ["_".join(sp.split()[:2]) for sp in dnm_df.Species]
    alpha_ped = {sp:[df["Pat DNMs"].sum(),df["Mat DNMs"].sum()] for sp,df in dnm_df.groupby("Species_simple")}
    dnm_ci = {}
    for sp in alpha_ped:
        paternal, maternal = alpha_ped[sp]
        alpha = paternal/maternal
        llimit, ulimit = statsmodels.stats.proportion.proportion_confint(paternal,
                                                                         paternal+maternal)
        dnm_ci[sp] = [llimit/(1-llimit),
                      alpha,
                      ulimit/(1-ulimit)
                     ]
    return dnm_ci

def model(mu_y, G_a, ratios, n_e, alpha, species, hue):
    '''Returns predicted alphas for range of sex-specific generation times'''
    alphas_r = []
    
    for r in ratios:
        G_f, G_m = ratioG(G_a,r)
        mu_f, mu_m = get_parental_age_effects_yrate(mu_y, G_f, G_m, n_e, alpha)
        a = predict_alpha(mu_f, mu_m, G_f, G_m, n_e)
        alphas_r.append(a)
        
    return pd.Series([n_e,(mu_y*G_a)-n_e,n_e/(mu_y*G_a),species,alphas_r[0], alphas_r[1], alphas_r[2], hue])

def parse_alpha_tsvs(experiment_string, alphas_directory):
    """Returns dataframe with alphas in all experiments found in comma-separated string"""
    merged_alphas_l = []
    for experiment in experiment_string.replace(" ", "").split(","):
        alpha_tsv = "{}/{}.LM.tsv".format(alphas_directory, experiment)
        alpha_df = pd.read_csv(alpha_tsv, sep="\t")
        alpha_df["experiment"] = experiment
        merged_alphas_l.append(alpha_df)
    merged_alphas_df = pd.concat(merged_alphas_l)
    return merged_alphas_df

def read_model_xy(g, directory):
    '''Return dictionary with model params and xy coordinates for given group'''
    model_params = pd.read_csv("{}/{}.model_params.tsv".format(directory, g),sep="\t")
    xy_data = pd.read_csv("{}/{}.xy_data.tsv".format(directory, g),sep="\t")
    return model_params, xy_data

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

def plot_pgls(data, g, experiment, pgls, ax, color, scatter, alpha, line, label_legend):
    '''Plots experiment vs. X(Z)/A alpha in group g to ax,
    also returns dict with lambda, r2 and pval'''

    # Get params into dict, then get pertinent stats
    param_df = data[g]["params"]
    params = list(param_df[param_df.experiment==experiment].transpose().to_dict().values())[0]
    slope = params["{}_slope".format(pgls)]
    intercept = params["{}_intercept".format(pgls)]
    lambda_val = params["ml_lambda"]
    r2 = params["{}_rsq".format(pgls)]
    p = params["{}_pval".format(pgls)]
    
    # Get xy data
    xy_df = data[g]["xy"]
    xy = xy_df[xy_df.experiment==experiment].reset_index(drop=True)
    
    # Plot scatter and line
    label = ""
    if label_legend:
        if label_legend==True:
            label = "{}:\n".format(g) + r"$\lambda={:.2f}, p={:.3f}, r^2={:.2f}$".format(lambda_val, p, r2)
        else:
            label = label_legend
    if scatter:
        ax.scatter(xy["xvar"], xy["yvar"], label=label, c=color, alpha=alpha, zorder=10)
    min_x = xy["xvar"].min()
    max_x = xy["xvar"].max()
    min_y = min_x*slope + intercept
    max_y = max_x*slope + intercept
    if line:
        ax.plot([min_x, max_x], [min_y, max_y], c=color)

    return {"lambda":lambda_val, "r2":r2, "pval":p}

def plot_sp_text(data, g, experiment, sp2common, selected_species, ax, colors):
    '''Annotates species names to ax'''

    # Get xy data
    xy_df = data[g]["xy"]
    xy = xy_df[xy_df.experiment==experiment].reset_index(drop=True)
    xy = xy[xy.Species.isin(selected_species)]
    texts = [ax.text(r.xvar, r.yvar, sp2common[r.Species], size=9, c=colors[i]) for i,r in xy.iterrows()]
    adjust_text(texts,arrowprops=dict(arrowstyle='-', color='gray'),expand_text=[2,2])

def model_alpha_vs_time(rate, n_e, Gs, alpha_p):
    '''Returns expected alpha for each timepoint in Gs, given
    a ratio of parental age effects alpha_p & n_e EE mutations'''

    alphas = []
    frac_p = alpha_p/(alpha_p+1)
    
    for g in Gs: 
        p_rate = (g*rate)
        p_rate = p_rate-n_e if p_rate-n_e>0 else 0
        frac = np.average([0.5, frac_p], weights=[n_e, p_rate])
        alpha = frac/(1-frac)
        alphas.append(alpha)
        #paternal = paternal_yearly_rate*g
        #maternal = maternal_yearly_rate*g
        #alphas.append((n_e+paternal)/(n_e+maternal))
    return alphas

def get_parental_age_effects_yrate(mu_y, G_f, G_m, n_e, alpha_p):
    '''Returns parental age effects given:
    (I)   Yearly autosomal substitution rate (mu_y)
    (II)  Sex-specific generation times (G_*)
    (III) Ratio of paternal age effects in (alpha_p)
    (IV)  Early embryonic mutations (n_e)'''
    
    n = (mu_y*(G_f + G_m))-(2*n_e)
    d = G_f + alpha_p*G_m
    mu_f = n/d
    mu_m = alpha_p*mu_f
    return mu_f, mu_m

def get_parental_age_effects(mu_g, G_f, G_m, n_e, alpha_p):
    '''Returns parental age effects given:
    (I)   Per generation rate (mu_g)
    (II)  Sex-specific generation times (G_*)
    (III) Ratio of paternal age effects in (alpha_p)
    (IV)  Early embryonic mutations (n_e)'''
    
    n = 2*mu_g - 2*n_e
    d = G_f + alpha_p*G_m
    mu_f = n/d
    mu_m = alpha_p*mu_f
    return mu_f, mu_m

def predict_alpha(mu_f, mu_m, G_f, G_m, n_e):
    '''Predicts alpha given:
    (I)   Parental age effects (mu_f, mu_m)
    (II)  Sex-specific generation times (G_*)
    (III) Early embryonic mutations (n_e)'''
    
    maternal = G_f*mu_f + n_e
    paternal = G_m*mu_m + n_e
    return paternal/maternal

def predict_yearly_rate(mu_ee, mu_f, mu_m, G_f, G_m):
    '''Predicts yearly mutation rate'''
    return (2*mu_ee + (mu_f*G_f) + (mu_m*G_m))/(G_f+G_m)

def ratioG(G_a, ratio):
    '''Returns sex-specific generation times given a male-to-female ratio'''
    G_f = (2*G_a)/(1+ratio)
    return G_f,G_f*ratio

def add_binomial_CIs(df):
    '''Adds binomial CIs to counts of Paternal vs. Maternal DNMs'''

    ci_df = pd.DataFrame()
    
    for i,r in df.iterrows():
        paternal, maternal = r.Paternal, r.Maternal
        ci = ssp.proportion_confint(paternal, paternal+maternal)
        estimates = [ci[0]/(1-ci[0]),paternal/maternal,ci[1]/(1-ci[1])]
        ci_df = pd.concat([ci_df, pd.DataFrame(estimates).transpose()])

    ci_df.columns = ["low","point","high"]
    ci_df = ci_df.reset_index(drop=True)
    return pd.concat([df, ci_df],axis=1)

def dodge_positions(df, var, order_var, order, margin):
    '''Slightly changes positions in xvar by margin units,
    keeping the central point in the same coordiante'''

    new_positions = []
    for i,r in df.iterrows():
        iorder = order.index(r[order_var])
        if iorder==0:
            new_positions.append(r[var]-margin)
        elif iorder==1:
            new_positions.append(r[var])
        elif iorder==2:
            new_positions.append(r[var]+margin)
    df["{}_dodge".format(var)] = new_positions
    return df

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    '''Truncates cmap for given range'''
    
    new_cmap = colors_mplt.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def extraploate_parental_dnms(df):
    '''Extrapolates number of parental DNMs assuming alpha_stage is constant'''
    
    p_dnms, m_dnms = 0,0
    for i,row in df.iterrows():
        p_frac = row.Paternal/(row.Maternal+row.Paternal)
        m_frac = 1 - p_frac
        p_dnms += p_frac*row.Total
        m_dnms += m_frac*row.Total
    return m_dnms, p_dnms

def fraction_ee(df):
    '''Returns fraction of DNMs in EE and periPGCS'''
    
    ee = df[df.Stage=="Early_embryonic"]["Total"].values[0]
    postpgcs = df[df.Stage=="PostPGCS"]["Total"].values[0]
    return ee/(ee+postpgcs)
