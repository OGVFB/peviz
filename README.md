# peviz
Shiny (Server) visualization of protein evolution (MSA)

# Function
Visualize [MSA](https://en.wikipedia.org/wiki/Multiple_sequence_alignment) across user specified proteins (via UniProt api) with accompany phylogenetic tree.

Optionally can display UniProt annotations

# Some (R or R adjacent) tools for MSA

  - https://github.com/plotly/react-msa-viewer
    - plotly based viewer
  - http://yulab-smu.top/ggmsa/index.html
    - "gg" like visualization
    - Can't use standard geoms? 
    - VERY slow and odd plot result (vertically squished) with large protein (Abca4; over 2000 AA)
  - http://www.bioinf.jku.at/software/msa/
    - MSA visualization done with LaTex
    - challenging to get in a Shiny / reactive environment as output is either pdf or Tex code
  - https://github.com/mhahsler/rBLAST/
    - R interface for blast
    - blast itself needs to be installed
  - https://github.com/vragh/seqvisr
  - https://github.com/GMOD/jbrowse-plugin-msaview 
    - jbrowse2 plugin....and jbrowseR package can display in Shiny
    - not certain if I can get msaview into Shiny via jbrowseR....
  - http://bioconductor.org/packages/release/bioc/html/DECIPHER.html
    - http://www2.decipher.codes/Gallery.html
    



