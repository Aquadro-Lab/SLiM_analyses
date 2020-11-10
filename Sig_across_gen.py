import pandas as pd
import os
from glob import iglob
from pathlib import Path
import re
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import sys
import argparse


#python sig_across_gen PATH, FIG_TITLE, WOL_INTERVAL, WOL_FIRST

PATH = sys.argv[1]
FIG_TITLE = sys.argv[2]
WOL_INTERVAL = sys.argv[3] #in the future, see if we can get this directly from the file name
WOL_FIRST = sys.argv[4] #same as prev comment

print("")
print("MKT result files are found in", str(PATH))
print("Figures will be called", str(FIG_TITLE))
print("The Wolbachia interval is", str(WOL_INTERVAL), "k generations")
print("Wolbachia infection occurs first in the cycle:", str(WOL_FIRST))
print("")

def scatter_all_gen(path_to_MKT_folder, fig_title): #function to create scatter plots across generations
    neutralrootdir = path_to_MKT_folder
    pathlist = Path(neutralrootdir).rglob('MKT_FWW_slim_bam_cds_*.out')

    MKT_df = pd.DataFrame(columns = ['Run','Gen', 'Dn', 'Ds', 'Pn', 'Ps', 'alpha05', 'pval05', 'DnDs', 'PnPs'])
    MKT_graph_df = pd.DataFrame(columns = ['Run','Gen', 'Dn', 'Ds', 'Pn', 'Ps', 'alpha05', 'pval05', 'DnDs', 'PnPs'])
    
    for path in pathlist:
        if os.stat(path).st_size != 0:
            data_temp = pd.read_table(path, header=None)
        
            MKT_headers = data_temp[data_temp[0].str.contains('Gen')]
        
            for i in range(len(MKT_headers)):
                MKT_headers_row = MKT_headers.index[i]
                MKT_headers_title = MKT_headers.values[i]
            
                gen_temp = re.search('._Gen(.*)k_.', MKT_headers_title[0])
                gen = int(gen_temp.group(1))
            
                run_temp = re.search('.k_(.*)', MKT_headers_title[0])
                run = run_temp.group(1)
            
                mktTable_text = data_temp.iloc[MKT_headers_row + 2]
            
                if mktTable_text[0] == "$mktTable":
                
                    gen_index = data_temp.iloc[MKT_headers_row]                 
                    neutral = data_temp.iloc[MKT_headers_row + 4]
                    selected = data_temp.iloc[MKT_headers_row + 5]
                    cutoff05 = data_temp.iloc[MKT_headers_row + 9]

                    neutral_split = neutral.str.split()
                    Ps = int(neutral_split[[0]][0][2])
                    Pn = int(neutral_split[[0]][0][3])

                    selected_split = selected.str.split()
                    Ds = int(selected_split[[0]][0][2])
                    Dn = int(selected_split[[0]][0][3])

                    cutoff05_split = cutoff05.str.split()
                    cutoff05_alpha = float(cutoff05_split[[0]][0][2])
                    cutoff05_pval = float(cutoff05_split[[0]][0][3])
                
                    DnDs = Dn/Ds
                    PnPs = Pn/Ps
            
                MKT_df = MKT_df.append(pd.Series([run, gen, Dn, Ds, Pn, Ps, cutoff05_alpha, cutoff05_pval, DnDs, PnPs], index=MKT_df.columns), ignore_index=True)
                MKT_graph_df = MKT_graph_df.append(pd.Series([run, gen, Dn, Ds, Pn, Ps, cutoff05_alpha, cutoff05_pval, DnDs, PnPs], index=MKT_graph_df.columns), ignore_index=True)
    
#    print(MKT_graph_df)

    MKT_graph_df.loc[MKT_df['pval05'] < 0.05, 'Sig'] = 'p<0.05'
    MKT_graph_df.loc[MKT_df['pval05'] >= 0.05, 'Sig'] = 'p>0.05'
    
    maxPnPs = max(MKT_graph_df['PnPs'])
    maxDnDs = max(MKT_graph_df['DnDs'])
    
    larger_axis = max(maxDnDs, maxPnPs)
    
    custom_palette = {"p<0.05":"darkorange", "p>0.05":"steelblue"}  

    sns.set(font_scale=2)

    group = "Gen"
    num_plots = len(MKT_df['Gen'].unique().tolist())
    total_rows = 4
    total_cols = int(num_plots/4 + (0 if num_plots % 4 == 0 else 1))

    fig, axs = plt.subplots(total_rows, total_cols, 
                            sharey=True, 
                            figsize=(7*total_cols, 7*total_rows),
                            constrained_layout=True)

    for (group, MKT_graph_df),ax in zip(MKT_graph_df.groupby(group),axs.flat):
        
#        sig_temp = MKT_graph_df.groupby(group)
        sig_count = (MKT_graph_df['Sig'] == 'p<0.05').sum()
        sig_text=AnchoredText("p<0.5 = " + str(sig_count), loc="upper right", frameon=False)
        
        
        ax.set_title(group)
        sns.scatterplot(data=MKT_graph_df, x="PnPs", y="DnDs", hue="Sig", palette=custom_palette, s=200, ax=ax, legend=False, alpha=0.6)
        ax.set_xlim(0,maxPnPs)
        ax.plot([0, larger_axis], [0, larger_axis], linewidth=1, color="darkgray")
        ax.add_artist(sig_text)
    
    fig.suptitle(fig_title, fontweight='bold')
    
#    pd.set_option('display.max_rows', 10) #this is for viewing the dataframe that is commented out below
    
#    return(MKT_df.sort_values(by='Gen')) #to view the dataframe
#    return(fig)
    fig_file_name = fig_title + "_scatter.png"
    fig_file_name = str(fig_file_name)
    plt.savefig(fig_file_name, format='png')
    plt.close()


def sig_across_gens_backcol(path_to_MKT_folder, fig_title, wol_interval, wolfirst): #function to create linegraph of significant counts across generations
    neutralrootdir = path_to_MKT_folder
    pathlist = Path(neutralrootdir).rglob('MKT_FWW_slim_bam_cds_*.out')

    MKT_df = pd.DataFrame(columns = ['Run','Gen', 'Dn', 'Ds', 'Pn', 'Ps', 'alpha05', 'pval05', 'DnDs', 'PnPs'])
    MKT_graph_df = pd.DataFrame(columns = ['Run','Gen', 'Dn', 'Ds', 'Pn', 'Ps', 'alpha05', 'pval05', 'DnDs', 'PnPs'])
    
    for path in pathlist:
        if os.stat(path).st_size != 0:
            data_temp = pd.read_table(path, header=None)
        
            MKT_headers = data_temp[data_temp[0].str.contains('Gen')]
        
            for i in range(len(MKT_headers)):
                MKT_headers_row = MKT_headers.index[i]
                MKT_headers_title = MKT_headers.values[i]
            
                gen_temp = re.search('._Gen(.*)k_.', MKT_headers_title[0])
                gen = int(gen_temp.group(1))
            
                run_temp = re.search('.k_(.*)', MKT_headers_title[0])
                run = run_temp.group(1)
            
                mktTable_text = data_temp.iloc[MKT_headers_row + 2]
            
                if mktTable_text[0] == "$mktTable":
                
                    gen_index = data_temp.iloc[MKT_headers_row]                 
                    neutral = data_temp.iloc[MKT_headers_row + 4]
                    selected = data_temp.iloc[MKT_headers_row + 5]
                    cutoff05 = data_temp.iloc[MKT_headers_row + 9]

                    neutral_split = neutral.str.split()
                    Ps = int(neutral_split[[0]][0][2])
                    Pn = int(neutral_split[[0]][0][3])

                    selected_split = selected.str.split()
                    Ds = int(selected_split[[0]][0][2])
                    Dn = int(selected_split[[0]][0][3])

                    cutoff05_split = cutoff05.str.split()
                    cutoff05_alpha = float(cutoff05_split[[0]][0][2])
                    cutoff05_pval = float(cutoff05_split[[0]][0][3])
                
                    DnDs = Dn/Ds
                    PnPs = Pn/Ps
            
                MKT_df = MKT_df.append(pd.Series([run, gen, Dn, Ds, Pn, Ps, cutoff05_alpha, cutoff05_pval, DnDs, PnPs], index=MKT_df.columns), ignore_index=True)
                MKT_graph_df = MKT_graph_df.append(pd.Series([run, gen, Dn, Ds, Pn, Ps, cutoff05_alpha, cutoff05_pval, DnDs, PnPs], index=MKT_graph_df.columns), ignore_index=True)
    
    MKT_graph_df.loc[MKT_df['pval05'] < 0.05, 'Sig'] = 'p<0.05'
    MKT_graph_df.loc[MKT_df['pval05'] >= 0.05, 'Sig'] = 'p>0.05'
    
    Gen_list = sorted(MKT_df['Gen'].unique().tolist())

    sig_only_df = MKT_graph_df[MKT_graph_df.Sig == 'p<0.05']    

    sig_only_df2 = sig_only_df.loc[sig_only_df['DnDs'] > sig_only_df['PnPs'], 'DnDs>PnPs'] = 'Yes'
    sig_only_df2 = sig_only_df.loc[sig_only_df['DnDs'] < sig_only_df['PnPs'], 'DnDs>PnPs'] = 'No'  
    sig_graph_df2 = sig_only_df.groupby(['Gen', 'DnDs>PnPs']).Gen.agg('count').to_frame('# p<0.05').reset_index()
    
    empty_rows = []
    empty_rows_cols = ['Gen', 'DnDs>PnPs', "# p<0.05"]
    for k in Gen_list:
        if ((sig_graph_df2['Gen'] == k) & (sig_graph_df2['DnDs>PnPs'] == 'Yes')).any() == False:
            empty_rows.append([k, 'Yes', 0])

        if ((sig_graph_df2['Gen'] == k) & (sig_graph_df2['DnDs>PnPs'] == 'No')).any() == False:
            empty_rows.append([k, 'No', 0])
    
    empty_rows_df=pd.DataFrame(empty_rows, columns=empty_rows_cols)
    
    sig_graph_df3 = sig_graph_df2.append(empty_rows_df)
    sig_graph_df4 = sig_graph_df3[sig_graph_df3['DnDs>PnPs'] == "Yes"]
    
#    print(sig_graph_df4)
    
    sns.set(rc={'figure.figsize':(12,3)})

#    linegraph = sns.lineplot(x="Gen", y="# p<0.05", hue="DnDs>PnPs", marker="o", data=sig_graph_df3)
    linegraph = sns.lineplot(x="Gen", y="# p<0.05", marker="o", data=sig_graph_df4)
    linegraph.set_title(fig_title, fontweight='bold')
    linegraph.axhline(1.25)
    plt.legend(fontsize='x-small', title_fontsize='6')
    
    def drange(start, stop, step):
        r=start
        while r< stop:
            yield r
            r += step
            
    if wolfirst == "True":
        backgroundcols = drange(0, max(Gen_list), 2*float(wol_interval))
    else:
        backgroundcols = drange(0+float(wol_interval), max(Gen_list), 2*float(wol_interval))
        
    for c in backgroundcols:
#        print(c)
        plt.axvspan(c, c+float(wol_interval), facecolor='y', alpha=0.2, zorder=-100)
    
#    return(linegraph)
    fig_file_name = fig_title + "_line.png"
    fig_file_name = str(fig_file_name)
    plt.savefig(fig_file_name, format='png')
    plt.close()


#scatter_all_gen('/mnt/c/Users/miwaw/Documents/Lab_Work/bam_seq_analysis/from_Runxi/10.22.20/bam_div_vs_con_conserv/slim_bam_cds_divAA_s1e-03_phase12pt5k/MKT_FWW_wolplus', "Conflict_s0.001_12.5k_wolFirst")
#sig_across_gens_backcol('/mnt/c/Users/miwaw/Documents/Lab_Work/bam_seq_analysis/from_Runxi/10.22.20/bam_div_vs_con_conserv/slim_bam_cds_divAA_s1e-03_phase12pt5k/MKT_FWW_wolplus', "Conflict_s0.001_12.5k_wolFirst", 12.5, True)

scatter_all_gen(PATH, FIG_TITLE)
sig_across_gens_backcol(PATH, FIG_TITLE, WOL_INTERVAL, str(WOL_FIRST))