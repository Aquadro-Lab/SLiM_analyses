# SLiM_analyses
Scripts to analyze output data from nucleotide based SLiM pipeline that models different modes of selection.

## Sig_across_gens.py
Takes in a folder with all the MKT_FWW text output files for a specific parameter set and outputs two graphs that display the number of significant MKT results across generation timepoints.  
  
**Use is:** Sig_across_gens.py path_to_folder_with_MKT_files title_for_figures wolbachia_infection_length wolbachia_first.  
**E.g.** Sig_across_gens.py /home/bam_div_vs_con_conserv/slim_bam_cds_divAA_s1e-03_phase12pt5k/MKT_FWW_wolplus Conflict_s0.001_12.5k_wolFirst 12.5 True  
**This would output the example files that are in this folder:** Conflict_s0.001_12.5k_wolFirst_line.png and Conflict_s0.001_12.5k_wolFirst_scatter.png.
  
**Ideas to improve:**  
1. Make it so the you only have to give the path name and the script can generate the figure title and identify the wolbachia_infection_length and Wolbachia_first parameters.  
2. Add a parameter so you can specify what the axes are in the line graph to make them easy to compare across different parameter setups. 
