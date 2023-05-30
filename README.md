## Working branch for web development

## Aim
**To develop a web framework written in python (Flask) for interactive plotting and visualisation of isoform dataset.**

- Gene level  
  -  Box-plot of the gene expression data grouped by developmental stage (pre- and post-natal) and sex
  -  User can select gene of choice
- Transcript level  
  -  Bar-plot of the number of FSM, ISM, NIC, NNC isoforms per gene 
  -  Tracks of the isoform (user can select isoform to show isoform expression)
  -  Box-plot of the isoform expression (also grouped by developmental stage and sex)

## Metadata
1. demux.csv = isoform expression (full-length read counts)
    - first column = list of isoforms
    - subsequent columns = sample IDs
    - such that each row is an isoform with the counts for each sample 

2. manifest.csv = SQANTI annotation of the isoforms 
    - further details for each isoform; for the purpose of this project, only refer to the `isoform`, `associated_gene` and `structural_category` column

3. phenotype.csv = phenotype data of the samples
    - sample IDs in phenotype.csv should match the sample IDs in demux.csv

## App
The web is written in Python3 and html using Flask.   
To run: `python main.py`

### Python-related libraries 
- Flask
- Pandas
- Plotly
