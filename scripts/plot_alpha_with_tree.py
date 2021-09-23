#!/usr/bin/env python

### Modules ###
import pandas as pd
from ete3 import Tree
import toyplot.pdf
import toytree
import sys
import statsmodels.stats.proportion
from optparse import OptionParser
import numpy as np

### Functions ###
def parse_alpha_tsvs(experiment_string, alphas_directory):
    # Returns dataframe with alphas in all experiments found in comma-separated string
    merged_alphas_l = []
    for experiment in experiment_string.replace(" ", "").split(","):
        alpha_tsv = "{}/{}.LM.tsv".format(alphas_directory, experiment)
        alpha_df = pd.read_csv(alpha_tsv, sep="\t")
        alpha_df["experiment"] = experiment
        merged_alphas_l.append(alpha_df)
    merged_alphas_df = pd.concat(merged_alphas_l)
    return merged_alphas_df

def discard_species(alphas, discard_file):
    # Returns dataframe without the species & experiments specified in discard_file
    discard_df = pd.read_csv(discard_file, sep="\t")
    ## Stderr
    for i,r in discard_df.iterrows():
        if r.species in alphas["species"].values:
            sys.stderr.write("Discarded {} from {}\n".format(r.species, r.experiment))
    discard_list = list(discard_df["species"] + "-" + discard_df["experiment"].values)
    alphas["experiment_species"] = alphas["species"] + "-" + alphas["experiment"]
    clean_alphas = alphas[~alphas.experiment_species.isin(discard_list)]
    clean_alphas.index = range(len(clean_alphas))
    return clean_alphas

def remove_bowels(w):
    # Remove lowercase bowels from word
    return "".join([l for l in w if l not in ["a","e","i","o","u"]])

def shorten(s, maxl=1000):
    # Shorten species label if too long, replacing capitalized words with starting letter + "."
    if len(s)<maxl:
        return s
    else:
        news = []
        olds = s.split()
        for i,w in enumerate(olds):
            if w[0].isupper():
                news.append(remove_bowels(w) + ".")
            currents = " ".join(news + olds[i+1:])
            print(currents, len(currents))
            if len(currents)<maxl:
                return(currents)
    return(currents)

def add_label_padding(sp_list, padding=4):
    # Add n white spaces (padding var) to the longest label, avoiding cropping by right-side of plot
    maxl = max([len(sp) for sp in sp_list])
    padchar = " "*padding
    return([sp+padchar if len(sp)==maxl else sp for sp in sp_list])

def miyata(xar):
    # Alpha estimate from X/A subsitution rate ratio
    return (4-(3*xar))/((3*xar)-2)

def reverse_miyata(alpha):
    # Given alpha, return X/A ratio
    return (2*alpha+4)/float(3*alpha+3)

def dodge_alphas(df):
    # Modifies the y axis to accomodate several observations per species
    if len(df)==1:
        df["position_tree"] = df["original_position_tree"].values
    else:
        df["position_tree"] = dodge_leaf_positions(df["original_position_tree"], 0.25)
    return df

def alpha_from_pedigrees(dnm_file):
    # Reads dnms from file and computes alpha and binomial CIs
    dnm_df = pd.read_csv(dnm_file, sep="\t")
    alpha_ped = {"_".join(r.Species.split()[:2]): [r["Pat DNMs"],r["Mat DNMs"]]
                 for i,r in dnm_df.iterrows()
                }
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

def style_tip_labels(maxl, labels):
    # New lines if tip labels are too long
    new_labels = []
    for sp in labels:
        if len(sp)>=maxl:
            splitsp = sp.split()
            sp = " ".join(splitsp[:2]) + "\n" + " ".join(splitsp[2:])
        new_labels.append(sp)
    return new_labels

def reduce_latin_labels(sp_list):
    # Shortens species names (eg. Homo_sapiens -> H.sapiens)
    return ["{}.{}".format(l.split("_")[0][0],l.split("_")[1]) for l in sp_list]

def convert_latin_to_common(sp_list, latin2common_df):
    # Converts latin names to common names
    l2c = latin2common_df.set_index("Species").to_dict()["CommonName"]
    return [l2c[sp] if sp in l2c else sp for sp in sp_list]

def dodge_leaf_positions(leaves, y_spacing):
    # Returns slightly modified leaf position (y axis),
    # so data from different experiments does not overlap
    n = len(leaves)
    dodge_space  = np.linspace(-y_spacing, y_spacing, n)
    return [original_y-extra_y for extra_y, original_y in zip(dodge_space, leaves)]

def emit_limit_stderr(observation, limit, sp, experiment):
    # Returns limit instead of observation and emits to stderr
    sys.stderr.write("Upper alpha limit for {} in {} has been reduced from {:.2f} to {}\n".format(sp,
                                                                                                  experiment,
                                                                                                  observation,
                                                                                                  limit
                                                                                                 )
                    )
    return limit

def replace_column_values(l, replace_dict):
    # Returns list of values with strings replaced according to dictionary
    return [element if element not in replace_dict else replace_dict[element] for element in l]

def plot_to_terminal(alphas, species):
    # Plots alpha point estimates to terminal
    sys.stdout.write("\n\nPlotting alpha point estimates to terminal:\n")
    fig = tpl.figure()
    fig.barh([alphas[sp] for sp in species], species, force_ascii=True)
    fig.show()
    sys.stdout.write("\n")
    
def restrict_alphas(alphas, upper, lower):
    # Restricts alphas to a given interval, assigns max/min 
    alphas["alpha"] = [upper if a>=upper else lower if a<=lower else a for a in alphas.alpha]
    alphas["lwr"] = [lower if l<lower else upper if l>upper else l for l in alphas.lwr]
    alphas["lwr_genratio"] = [lower if l<lower else upper if l>upper else l for l in alphas.lwr_genratio]
    alphas["upr"] = [upper if u>upper else lower if u<lower else u for u in alphas.upr]
    alphas["upr_genratio"] = [lower if l<lower else upper if l>upper else l for l in alphas.upr_genratio]
    return alphas

def transform_bounds(hb, vb):
    # Transform fraction bounds into toyplot lists
    hbounds = hb.split(",")
    vbounds = vb.split(",")
    ax0 = tuple(hbounds[:2] + vbounds)
    ax1 = tuple(hbounds[2:] + vbounds)
    return ax0, ax1

def plot_tree_alphas(alphas, nwk, labels, outpdf, width, height, legend, alpha_upper_limit, alpha_lower_limit, dnms, hue, horitzontal_bounds, vertical_bounds, palette):
    # Plots phylogenetic tree and aligned alphas
    
    ## If height of figure is not enough for xlabel, increase
    height = height if height>75 else height+50
    xlabel_margin = 100 - (50/height*100)
    vert_limits = list(map(int,vertical_bounds.replace("%","").split(",")))
    if vert_limits[1] > xlabel_margin:
        vertical_bounds = "{}%,{}%".format(vert_limits[0],xlabel_margin)
    
    ## Define 2 plots, tree & scatter with male bias
    canvas = toyplot.Canvas(width=width, height=height)
    ax0bounds, ax1bounds = transform_bounds(horitzontal_bounds, vertical_bounds)
    ax0 = canvas.cartesian(bounds=(ax0bounds))
    ax1 = canvas.cartesian(bounds=(ax1bounds))

    ## Plot tree
    style = {"edge_style":{"stroke-width": 1}, "tip_labels_style":{"font-size": "11px", "-toyplot-anchor-shift": "3px"}}
    nwk.draw(tip_labels_align=True, axes=ax0, tip_labels = labels, **style)
    ax0.show = False

    ## Define colors for several experiments
    if palette=="Paired":
        col_list = [c for i,c in enumerate(toyplot.color.brewer.palette(palette)) if i!=10] # Skip yellow, can't be seen in white background
        palette_p = {var:col_list[i] for i,var in enumerate(alphas[hue].unique())}
    if palette!="Paired":
        with open(palette, "r") as palette_file:
            palette_p = {line.split()[0]:line.split()[1] for line in palette_file if "Color" not in line}

    seen_line_to_leaves = []
    point_alphas = {}

    ## Restrict point alphas to interval
    alphas = restrict_alphas(alphas, alpha_upper_limit, alpha_lower_limit)

    # Plot alphas
    legends = {}
    for var,df in alphas.groupby(hue):
        legends[var] = ax1.scatterplot(df.alpha, df.position_tree,
                                              color=palette_p[var.split(".")[0]],
                                              size=3
                                             )
        for i,r in df.iterrows():
            leaf = r.position_tree
            of_leaf = r.original_position_tree
            lower = r.lwr
            upper = r.upr
            gen_lower = r.lwr_genratio
            gen_upper = r.upr_genratio

            # CIs
            ax1.plot([lower, upper], [leaf, leaf],
                     color=palette_p[var.split(".")[0]],
                     style={"stroke-width": 2},
                     opacity = 1
                    )

            # Uncertainty due to differences in sex-specific generation time
            ax1.plot([gen_lower, gen_upper], [leaf, leaf],
                     color=palette_p[var.split(".")[0]],
                     style={"stroke-width": 2},
                     opacity = 0.6
                    )

            # Line connecting to tree tips, only done once
            if r.species not in seen_line_to_leaves:
                ax1.plot([0, gen_lower], [of_leaf, of_leaf],
                         color="gray",
                         style={"stroke-width": 0.4, "stroke-dasharray":"2, 2"},
                         opacity = 0.7
                        )
                # Also take advantage that this is done only once per species to plot DNM CIs, if available
                if r.species in dnms:
                    lower_dnm, alpha, upper_dnm = dnms[r.species]
                    if upper_dnm>=alpha_upper_limit:
                        upper_dnm = emit_limit_stderr(upper_dnm, alpha_upper_limit, r.species, "DNMs")
                    yrectmargin = 0.3
                    ax1.rectangle(lower_dnm,upper_dnm,of_leaf-yrectmargin,of_leaf+yrectmargin,opacity=0.15,color="gray")
                    ax1.rectangle(alpha-0.02,alpha+0.02,of_leaf-yrectmargin,of_leaf+yrectmargin,opacity=0.30,color="red")
                point_alphas[labels[of_leaf]] = round(r.alpha, 1)
            seen_line_to_leaves.append(r.species)



    # Legend
    if legend:
        canvas.legend(list(legends.items()),
                      corner=("right", 50, "20%", "20%"),
                     )

    # Some ugly but necessary steps to align the tree and scatter
    extra = len(labels)*0.044*(250*(1/height))
    ax1.scatterplot(0,-extra, color="white")
    ax1.scatterplot(0, nwk.ntips-1+extra, color="white")
    print(alpha_upper_limit)
    ax1.scatterplot(alpha_upper_limit,0, color="white")

    # Line in alpha = 1
    ax1.rectangle(1-0.01,1+0.01,-extra,nwk.ntips-1+extra,opacity=0.3,color="gray")

    ax1.y.show = False
    ax1.x.ticks.show = True
    ax1.x.label.text = "\u03B1"

    # Pdf if specified
    if outpdf!="NaN":
        toyplot.pdf.render(canvas, outpdf)

    return point_alphas

### Options ###
usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-x", "--experiments", action="store", type="str", dest="experiments",help="Comma separated experiment codes")
parser.add_option("-d", "--discard_file", action="store", type="str", dest="discard_file",help="File with 2 columns with species to discard (Species\tExperiment)")
parser.add_option("-o", "--output_pdf", action="store", type="str", dest="pdf_file",help="Path to pdf output")
parser.add_option("-t", "--terminal_plot", action="store_true", default=False, dest="terminal_plot",help="Plot rudimentary plot to terminal with alpha point estimates")
parser.add_option("-c", "--color_by", action="store", type="str", dest="colorby",help="Which variable to use for colouring (experiment/mut_type)")
parser.add_option("-i", "--include", action="store", type="str", dest="include",help="Which mutation classes to include")
parser.add_option("-g", "--gen_interval", action="store", type="str", dest="gen_interval",help="Sex-specific generation time ratios to add uncertainty to the CIs")
parser.add_option("-n", "--nwk", action="store", type="str", dest="nwk",help="Phylogeny to get the topology from")
parser.add_option("-y", "--height_factor", action="store", type="float", dest="height_factor",help="Plot height")
parser.add_option("-l", "--plot_legend", action="store_true", dest="plot_legend",default=False,help="Plot legend?")
parser.add_option("-p", "--palette", action="store", type="str", dest="palette",help="Color palette")
#parser.add_option("-s", "--shorten_species_labels", action="store", type="int", dest="max_species_string_length",help="Maximum numbers of characters in any species label, shorten if above")

(options, args) = parser.parse_args()

### Check options ### 
if not options.experiments:
    sys.exit("Provide a comma separated list of experiment codes")
if not options.nwk:
    sys.exit("Provide a nwk tree to draw the topology from")
if options.terminal_plot:
    import termplotlib as tpl
pdf_file = "NaN"
if options.pdf_file:
    pdf_file = options.pdf_file
height_factor = 15
if options.height_factor:
    height_factor = options.height_factor

### Process data ###
experiments = options.experiments
alpha_dir = "alphas"
zoonomia_nwk = options.nwk
taxonomy_file = "../data/latin2common.txt"
dnm_file = "../data/dnm_est_noDups.tsv"

# Read input and remove species if outgroups
merged_alphas = parse_alpha_tsvs(experiments,
                                 alpha_dir)
if options.discard_file:
    merged_alphas = discard_species(merged_alphas, options.discard_file)
all_species = list(merged_alphas["species"].unique())

# Read phylogenies and remove species (e.g. outgroups)
nwktree = toytree.tree(zoonomia_nwk, tree_format=1)
species_to_remove = [sp for sp in nwktree.get_tip_labels() if sp not in all_species]
nwktree = nwktree.drop_tips(species_to_remove)
species_tips = nwktree.get_tip_labels()
leaf_positions_in_tree = {sp:i for i,sp in enumerate(nwktree.get_tip_labels())}
merged_alphas["original_position_tree"] = [leaf_positions_in_tree[sp] for sp in merged_alphas["species"]]

# Replace latin names for common names, add padd space so labels are not cropped
common_names_df = pd.read_csv(taxonomy_file, sep="\t")[["Species", "Common_names"]]
common_names = common_names_df.set_index("Species").to_dict()["Common_names"]
with open("../data/Species_to_chromosomes.txt","r") as chrom_level_file:
    chrom_level = [line.split()[0] for line in chrom_level_file]
species_labels = [common_names[sp] if sp not in chrom_level else "{}*".format(common_names[sp]) for sp in species_tips]
species_labels = add_label_padding(species_labels, padding=25)

# Restrict to mutation types that have been selected, and change names 
merged_alphas.index = range(len(merged_alphas))
if options.include:
    include = options.include.split(",")
    merged_alphas = merged_alphas[merged_alphas.mut_type.isin(include)]
replace_mutypes = {"exptotsub":"All"}
merged_alphas["mut_type"] = replace_column_values(merged_alphas["mut_type"], replace_mutypes)
merged_alphas.index = range(len(merged_alphas))

# Add uncertainty due to different sex-specific generation times
if not options.gen_interval:
    sys.exit("Please provide a 'minimum, max' generation time ratio. Use 1,1 to not add uncertainty to the CIs")
ratio_interval = list(map(float, options.gen_interval.split(",")))
merged_alphas["lwr_genratio"] = merged_alphas["alpha_lwr"].apply(lambda x: miyata(reverse_miyata(x)/reverse_miyata(ratio_interval[1])))
merged_alphas["upr_genratio"] = merged_alphas["alpha_upr"].apply(lambda x: miyata(reverse_miyata(x)/reverse_miyata(ratio_interval[0])))

# Do some changes to the dataframe so it can be represented in the tree
merged_alphas = merged_alphas.groupby("species").apply(lambda x: dodge_alphas(x))
dnm_ci = alpha_from_pedigrees(dnm_file)
merged_alphas.columns = ["alpha","lwr","upr"] + list(merged_alphas.columns[3:])

# Plot with tree
point_alphas = plot_tree_alphas(alphas = merged_alphas,
                                nwk = nwktree,
                                labels = species_labels,
                                outpdf = pdf_file,
                                width = 700,
                                height = height_factor*len(species_labels),
                                legend = options.plot_legend,
                                alpha_upper_limit = 6,
                                alpha_lower_limit = 0,
                                dnms = dnm_ci,
                                hue=options.colorby,
                                horitzontal_bounds = "5%,50%,55%,95%" if options.plot_legend==False else "5%,42%,46%,75%",
                                vertical_bounds = "5%,95%",
                                palette = "Paired" if not options.palette else options.palette
                                )
# Plot to terminal?
if options.terminal_plot:
    plot_to_terminal(point_alphas, species_labels[::-1])
