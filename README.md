# Transcription_Factor_Prioritization
Uses ReMap ChIP-seq db, publicly available expression data, and expression/luciferase experiments to determine transcription factors that likely regulate mLair1, hLAIR1, and hLAIR2.

Input datasets: bedtools intersect of our luciferase regions with ReMap datasets (supplied for mouse and human)

Publicly available datasets used in this project can be downloaded from GEO and ReMap:
1.  ReMap database of transcription factor binding nonredundant peaks (remap.univ-amu.fr)
2.  Immgen Human PBMCs RNAseq (NCBI GEO dataset GSE227743)30
3.	Immgen Mouse ULI Systemwide RNA-seq profiles (NCBI GEO dataset GSE109125)
4.	RNAseq profiling of murine BMDM (NCBI GEO dataset GSE199128)31 
5.	RNAseq profiling of human monocytes and MDM (NCBI GEO dataset GSE147311) (PMID 33761354)
6.	Immgen Mouse Cytokines: Interferons Affymetrix Mouse Gene 1.0 ST Array (NCBI GEO dataset GSE75306) 
