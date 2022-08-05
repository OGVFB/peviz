# peviz
Visualization of protein evolution (MSA)

# Goal
For a *single* gene (protein) visualize [MSA](https://en.wikipedia.org/wiki/Multiple_sequence_alignment) across user specified species. The MSA is likely most commonly done by running web-based tooling from either [NCBI](https://www.ncbi.nlm.nih.gov/projects/msaviewer/) or [EBI](https://www.ebi.ac.uk/Tools/msa/).

I would like to also overlay known protein domains (via uniprot?) over the MSA to help with quick(er) interpretation of the substance of evolutionary differences.

A phylogenetic tree would also be useful (I think) to see the evolutionary relationships ordered by branches. 

# Wrinkle
This project involves some poorly/non annotated genomes/transcriptomes

# Some (R or R adjacent) tools

  - https://github.com/plotly/react-msa-viewer
    - plotly based viewer
  - http://yulab-smu.top/ggmsa/index.html
    - "gg" like visualization
    - Can't use standard geoms? 
  - http://www.bioinf.jku.at/software/msa/
    - MSA visualization done with Tex  

# Workflow?:

1. Input gene name and for non-"reference" genomes the GFF(s)
2. Find protein ID
3. Extract protein sequence across user specified (?) well annotated species
4. Find most likely ortholog in poorly/non annotated species (via gff?) and extract protein sequence
5. Run MSA
6. Extract known protein domains
7. Produce report:
  - MSA + protein domains
  - Phylogenetic tree


