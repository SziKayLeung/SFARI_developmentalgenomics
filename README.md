## Working branch for web development

## Aim
**To develop a web framework written in python (Mercury) for interactive plotting and visualisation of isoform dataset.** 

User can select gene (gene name) of choice:  
- Basic summary information
  - Number of known and novel transcripts
- Gene level  
  -  Scatter-plot (R ggplot) of the gene expression data across development (pre- and post-natal) and coloured by sex
  -  Differential gene expression results
- Transcript level  
  - Differential transcript expression results across development and between sex
  - Scatter plots (R ggplot) of selected transcript expression data across development (pre- and post-natal) and coloured by sex
  - Tracks of selected DTE isoform (user can select isoform to show isoform expression)

All files are deposited under `0_mercury` directory.
`0_flask` is redundant.

## Metadata

1. phenotype.csv = phenotype data of the samples
    - sample IDs in phenotype.csv should match the sample IDs in demux.csv

## App
The web is written in Python3 and html using Mercury.   
To run: 
```
# activate conda environement
module load Miniconda2
conda env create --file mercury_environment.yml
source activate mercury

# run mercury
mercury add .\LR_Resource.ipynb
mercury run LR_Resource.ipynb
```

### Python-related libraries 
- Pandas
- numpy
- plotnine
- scipy
