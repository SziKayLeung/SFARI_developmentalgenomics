#!/usr/bin/env python

import pandas as pd
import numpy as np
import plotnine
from scipy import misc
from plotnine import *
from IPython.display import display, HTML

def subset_file(file, num, pattern):
    header = []
    matching_lines = []
    with open(file, 'r') as file:
        for ele, line in enumerate(file):
            if ele == 0:
                #line = line.replace('\n','')
                header.append(line)
            if pattern == line.split(",")[num]:
                line = line.replace('\n','')
                matching_lines.append(line)
    results = pd.DataFrame([l.split(",") for l in matching_lines])
    if len(results) > 0:
        results.columns = [h.split(",") for h in header][0]
    return(results)

def subset_file_num(file, num, pattern):
    header = []
    matching_lines = []
    with open(file, 'r') as file:
        for ele, line in enumerate(file):
            if ele == 0:
                line = line.replace('\n','')
                header.append(line)
            else:                    
                if pattern == int(line.split(",")[num]):
                    line = line.replace('\n','')
                    matching_lines.append(line)
    results = pd.DataFrame([l.split(",") for l in matching_lines])
    if len(results) > 0:
        results.columns = [h.split(",") for h in header]
    return(results)

def plot_expression(dat):
    # Create a DataFrame 'dat' and convert 'group' to a categorical variable
    dat['group'] = pd.Categorical(dat['group'], categories=['Prenatal', 'Postnatal'], ordered=True)
    dat["age"] = pd.to_numeric(dat["age"])
    dat["normalised_counts"] = pd.to_numeric(dat["normalised_counts"])
    
    # Define age scaling functions
    def fetal_ages_scale(age):
        return np.interp(age, (0, 40), (0, 40))

    def child_ages_scale(age):
        return np.interp(age, (0, 40), (41, 60))

    def adult_ages_scale(age):
        return np.interp(age, (40, 100), (75, 100))
    


    # Apply the scaling functions to create 'age.rescale' column
    dat['age.rescale'] = np.nan
    dat['age.rescale'] = np.where(dat['group'] == 'Prenatal', dat['age'].apply(fetal_ages_scale), dat['age.rescale'])
    dat['age.rescale'] = np.where((dat['group'] == 'Postnatal') & (dat['age'] <= 40), dat['age'].apply(child_ages_scale), dat['age.rescale'])
    dat['age.rescale'] = np.where((dat['group'] == 'Postnatal') & (dat['age'] > 40), dat['age'].apply(adult_ages_scale), dat['age.rescale'])
    
    breaks = [dat['age.rescale'].iloc[i] for i in [0, 14, 30, 31, 32, 39, 19]]
    labels = [dat['age'].iloc[i] for i in [0, 14, 30]] + ['40/0'] + [dat['age'].iloc[i] for i in [32, 39, 19]]

    width, height = 15, 6 
    x_vline = 40

    # Create the plot
    gg = (
        ggplot(aes(x='age.rescale', y='normalised_counts'), data=dat) +
        geom_point(aes(fill = 'sex', colour = 'sex')) +
        labs(x=None, y='log10 normalized counts') +
        geom_smooth(method = "loess") +
        scale_x_continuous(breaks=breaks, labels=labels) +
        theme_classic() + 
        theme(figure_size=(width, height)) +
        geom_vline(xintercept=16,linetype="dashed")
    )

    x1 = dat['age.rescale'].min() + 4
    x2 = dat['age.rescale'].min() + 50
    y = dat['normalised_counts'].max() + 0.3

    # Add text annotations for "Pre-natal" and "Post-natal"
    gg = gg + annotate("text", x =x1, y = y, label = "Pre-natal") + annotate("text", x =x2, y = y, label = "Post-natal")

    return(gg)

def transcript_output(classFile, phenotype, inputDTE, inputDTEnorm, pattern):
    
    dat1 = classFile[classFile["isoform"] == pattern].set_index("isoform")
    dat2 = inputDTE[inputDTE["isoform"] == pattern].set_index("isoform")[["log2FoldChange","pvalue","padj"]]
    x = dat2.join(dat1, how='left') 
    xselect = x.loc[[pattern]]
    display(HTML("<div>" + xselect.style.render() + "</div>"))
        
    norm = subset_file(inputDTEnorm, 0, pattern).T.drop(index="isoform")
    norm.columns = ["normalised_counts"]
    norm = pd.merge(norm, phenotype, left_index=True, right_index=False, right_on = "sample")
    transcriptPlot = plot_expression(norm)
    print(transcriptPlot)
    
    return(transcriptPlot)
    
def extract_gtf(pattern, gtfInput, type):
    matching_lines = []
    with open(gtfInput, 'r') as file:
        for line in file:
            if type == "reference":            
                if pattern == line.split(",")[1]:
                    line = line.replace('\n','')
                    matching_lines.append(line)
            else:
                if pattern == line.split(",")[0]:
                    line = line.replace('\n','')
                    matching_lines.append(line)
    gtfExtract = pd.DataFrame([l.split(",") for l in matching_lines])
    if type == "reference":
         gtfExtract.columns = ["transcript_id","gene_id","seqnames","strand","start","end"] 
    else:
         gtfExtract.columns = ["transcript_id","gene_id","seqnames","strand","start","end"] 
    gtfExtract["start"] = [int(i) for i in gtfExtract["start"]]
    gtfExtract["end"] = [int(i) for i in gtfExtract["end"]]
    return(gtfExtract)

def get_type(transcript_id):
    if "ENST" in transcript_id:
        return 'Reference'
    else:
        return 'LR'
    
def plot_structure(gtf, refgtf, transcript, gene):
    gtfExtract = extract_gtf(transcript, gtf,"LR")
    rgtfExtract = extract_gtf(gene, refgtf,type="reference")
    merged = pd.concat([gtfExtract,rgtfExtract])
    merged['Type'] = merged['transcript_id'].apply(get_type)
    return(merged)
                                   

