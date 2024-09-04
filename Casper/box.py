import os;
from sys import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn
import matplotlib.colors as mcolors


fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))
def read_perc(tool, file_dir):
    rv_perc = []
    sh_perc = []
    
    for diff in diffs:
        rv = []
        sh = []
        INT_rv = []
        INT_sh = []
        #for line in open('new_'+tool+'_gene_perc_'+str(diff)+'.txt', "r"):
        for line in open(file_dir+'filtered_'+tool+'_count_perc_'+str(diff)+'.txt', "r"):
            if line.strip().split(" ")[0] == "ID":
                continue
            if line.strip().split("\t")[0] == "Shuffled Average:":
                break
            if len(line.strip().split("\t")) > 1:
                if line.strip().split("\t")[1] == "_rv":
                    rv.append(float(line.strip().split("\t")[2])*100)
                    if diff == 0:
                        INT_rv.append(float(line.strip().split("\t")[5])*100)
                elif line.strip().split("\t")[1] == "_s1" or "_s2" or "_s3":
                    sh.append(float(line.strip().split("\t")[2])*100)
                    if diff == 0:
                        INT_sh.append(float(line.strip().split("\t")[5])*100)
            else:
                continue
        rv_perc.append(rv)
        sh_perc.append(sh)
        # if diff == 0:
        #     rv_perc.append(INT_rv)
        #     sh_perc.append(INT_sh)
    return rv_perc,sh_perc


def draw_rv(rvs, buffer,color):
    legend_handles1 = []
    for i in range(0,len(tools)):
        rv_perc = rvs[i]
        pos = np.arange(len(rv_perc))+buffer[i]
        col = color[i]
        box = axes[0].boxplot(rv_perc, positions = pos,
                        boxprops=dict(facecolor=col, alpha = 0.5), widths = 0.13, vert=True, patch_artist=True , flierprops={'marker': '.'})
        # for i, data in enumerate(rv_perc):
        #     x = [pos[i] for t in range(len(data)) ]
        #     axes[0].scatter(x, data, color='black', marker='o', alpha=0.5)

        legend_handles1.append(box['boxes'][0])
        
        # Add the combined legend to the plot
    axes[0].legend(legend_handles1, tools, loc='lower right')

    axes[0].set_title('Reversed Complement')
    axes[0].set_ylabel('percentage of consistency %')
    axes[0].set_ylim([30,90])
    axes[0].set_xlabel('threshold')
    #axes[0].set_xticklabels(['exact match','integer','within 1%','within 10%'])
    num_ticks = 3
    xtick_pos = np.linspace(0, len(rv_perc) - 1, num_ticks)
    axes[0].set_xticks(xtick_pos)
    axes[0].set_xticklabels(['exact match','within 1%','within 10%'])
    for pos in xtick_pos[1:]:
        axes[0].axvline(x=pos - 0.5, linestyle='--', color='gray', alpha=0.7)

    xtick_pos = np.arange(len(rv_perc))
    axes[0].set_xticks(xtick_pos)
    

    
def draw_sh(shs,buffer,color):
    legend_handles1 = []
    for i in range(0,len(tools)):
        sh_perc = shs[i]
        pos = np.arange(len(sh_perc))+buffer[i]
        col = color[i]
        box = axes[1].boxplot(sh_perc, positions = pos,
                        boxprops=dict(facecolor=col, alpha = 0.5), widths = 0.13, vert=True, patch_artist=True, flierprops={'marker': '.'})
        # for i, data in enumerate(sh_perc):
        #     x = [pos[i] for t in range(len(data)) ]
        #     axes[1].scatter(x, data, color='black', marker='o', alpha=0.5)
        legend_handles1.append(box['boxes'][0])
    axes[1].legend(legend_handles1, tools, loc='upper left')
    axes[1].set_title('Shuffled')
    axes[1].set_ylabel('percentage of consistency %')
    axes[1].set_ylim([30,90])
    axes[1].set_xlabel('threshold')
    #axes[1].set_xticklabels(['exact match','integer','within 1%','within 10%'])
    num_ticks = 3
    xtick_pos = np.linspace(0, len(sh_perc) - 1, num_ticks)
    axes[1].set_xticks(xtick_pos)

    axes[1].set_xticklabels(['exact match','within 1%','within 10%'])
    for pos in xtick_pos[1:]:
        axes[1].axvline(x=pos - 0.5, linestyle='--', color='gray', alpha=0.7)

    xtick_pos = np.arange(len(sh_perc))
    axes[1].set_xticks(xtick_pos)
    

def load_and_plot_data(tools, file_dir):
    rvs = []
    shs = []

    # Adjust colors and buffers based on the number of tools
    num_tools = len(tools)
    # color = plt.cm.viridis(np.linspace(0, 1, num_tools))
    # color = plt.cm.tab10(np.linspace(0, 1, num_tools))
    # color = plt.cm.Set1(np.linspace(0, 1, num_tools))
    basic_colors = list(mcolors.TABLEAU_COLORS.values())
    color = basic_colors[:num_tools]
    buffer = np.linspace(-0.3, 0.3, num_tools)

    for i, tool in enumerate(tools):
        rv, sh = read_perc(tool, file_dir)
        rvs.append(rv)
        shs.append(sh)

    # buffer = buffer + np.arange(num_tools) * 0.1  # Adjust the buffer for better visualization

    draw_rv(rvs, buffer, color)
    draw_sh(shs, buffer, color)
    print(rvs)
    print(shs)
    plt.suptitle("Percentage of consistency between the original sample and replicates on gene count", wrap=True)
    plt.savefig("boxplot_nointeger_synthetic.png", dpi=199)
    plt.show()

if __name__ == "__main__":
    diffs = [0,0.01,0.1]
    tools = ['kallisto']#, 'Salmon', 'HTseq', 'IsoEm2', 'featurecounts', 'Rsem', 'STAR']  # Add or remove tools as needed
    file_dir = "/Users/fatemehmohebbi/Downloads/rna_seq_results/consistency_files_synthetic/"

    load_and_plot_data(tools, file_dir)


"""if __name__ == "__main__":
    diffs = [0,0.01,0.1]
    tools = ['kallisto','salmon', 'htseq', 'IsoEm2']#, 'IsoEm2', 'featurecounts']
        
    file_dir = "/Users/fatemehmohebbi/Downloads/rna_seq_results/consistency_files/"
    ka_rv,ka_sh = read_perc('kallisto', file_dir)
    sa_rv,sa_sh = read_perc('salmon', file_dir)
    ht_rv,ht_sh = read_perc('htseq', file_dir)
    is_rv,is_sh = read_perc('IsoEm2', file_dir)
    
    rvs = [ka_rv, sa_rv, ht_rv, is_rv]
    shs = [ka_sh, sa_sh, ht_sh, is_sh]
    buffer = [-0.3, -0.1, 0.1, 0.3]
    color = ['C0','salmon', 'm', 'r']
        
    draw_rv(rvs, buffer, color)
    draw_sh(shs, buffer, color)
    plt.suptitle("percentage of consistency between the original sample and replicates on gene count",
             wrap = True)
    #plt.savefig("boxplot.png", dpi=199)
    plt.savefig(file_dir+"boxplot_nointeger.png", dpi=199)"""


