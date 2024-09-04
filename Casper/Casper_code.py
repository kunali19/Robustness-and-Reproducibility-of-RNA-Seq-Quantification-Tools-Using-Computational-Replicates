from sys import *
import os;
os.environ["OMP_NUM_THREADS"] = "1"
import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np

#IDs = {"ERR188218", "ERR188434", "ERR188444"}
IDs = {"ERR188218", "ERR188434", "ERR188444"}
#A = {"SRR896663", "SRR896679", "SRR896695", "SRR896711", "SRR896727"}

# Address = "/Users/fatemehmohebbi/Downloads/rna_seq_results/salmon_output/salmon_"
# Ending = "_output/quant.sf"
# Address = "/Users/fatemehmohebbi/Downloads/rna_seq_results/kallisto_output/kallisto_"
# Ending = "_output/abundance.tsv"
# Address = "/Users/fatemehmohebbi/Downloads/rna_seq_results/IsoEM2_output/"
# Ending ="_bowtie2/Genes/gene_tpm_estimates"
# Address="/Users/fatemehmohebbi/Downloads/rna_seq_results/htseq_synthetic_data/"
# Ending="_output_htseq.txt"
# Address = "/Users/fatemehmohebbi/Downloads/rna_seq_results/rsem_synthetic_results/"
# output_dir = "/Users/fatemehmohebbi/Downloads/rna_seq_results/rsem_synthetic_results/tpm_resluts/"
# Ending = "_rsem.genes.resultstpm.txt"
# key = "Rsem"

Address = "/project/fmohebbi_1178/RNA-seq_project/CASPER/results/"
output_dir = "/project/fmohebbi_1178/RNA-seq_project/CASPER/final_output/"
Ending = ".txt"
key = "Rsem"

tpm_pos = -1

skip = "gene_id" #salmon
# skip = "target_id" #kallisto
# key = "IsoEM2"
# key = "htseq"

diffs = [0,0.01,0.1]
for diff in diffs:
    AVGDiagonal = 0
    AVGDiagonal_s1 = 0
    AVGDiagonal_s2 = 0
    AVGDiagonal_s3 = 0
    AVGXaxis = 0
    AVGYaxis = 0
    RVDia = 0
    RVX = 0
    RVY = 0
    diagonals = []
    rvs = []
    if os.path.exists(output_dir+'filtered_'+key+'_TPM_perc_'+str(diff)+'.txt'):
        os.remove(output_dir+'filtered_'+key+'_TPM_perc_'+str(diff)+'.txt')
    with open(output_dir+'filtered_'+key+'_TPM_perc_'+str(diff)+'.txt', 'w') as fi:
        if diff == 0:
            fi.write("ID Version DiagonalPercent XAxisPercent YAxisPercent Integer\n")
        else:
            fi.write("ID Version DiagonalPercent XAxisPercent YAxisPercent\n")
        for ID in IDs:
            OrginalResults = Address+ID+"_g"+Ending #.sf file
            versions = {"_rv","_s1","_s2","_s3"}
            OGLines = {}
            for line in open(OrginalResults, "r"):
                GeneName = line.strip().split("\t")[0]
                TPM = line.strip().split("\t")[tpm_pos]
                if GeneName == skip: #skip the first line
                   continue
                OGLines[GeneName] = TPM
            #versions = {"_rv"}
            for v in versions:
                TranscriptsOnTheDiagonal = 0
                Int_Diagonal = 0
                TotalTranscripts = 0
                ZeroInCase1 = 0
                ZeroInCase2 = 0
                MDLines = {}

                for line in open(Address+ID+v+Ending, "r"):
                    GeneName = line.strip().split("\t")[0]
                    TPM = line.strip().split("\t")[tpm_pos]
                    if GeneName == skip:
                      continue
                    MDLines[GeneName] = TPM

                SharedOGTPMs = []  # Lines in the Salmon results that are also in the Salmon results
                SharedModifiedTPMs = []  # Lines in the Salmon results that are also in the Salmon results
                TPM1 = 0
                TPM2 = 0
                for GeneName in OGLines:
                    if GeneName in MDLines:
                        TPM1 = float(OGLines[GeneName])
                        TPM2 = float(MDLines[GeneName])

                        if TPM1 * (1+diff) >= TPM2 >= TPM1 * (1-diff):
                            TranscriptsOnTheDiagonal += 1

                        if TPM1 <= 0.0 and TPM2 > 0.1:
                            ZeroInCase1 += 1

                        if TPM2 <= 0.0 and TPM1 > 0.1:
                            ZeroInCase2 += 1
                        if diff == 0:
                            if int(TPM1) == int(TPM2):
                                Int_Diagonal += 1
                        SharedOGTPMs.append(math.log(TPM1+1, 10))
                        SharedModifiedTPMs.append(math.log(TPM2+1, 10))
                        TotalTranscripts += 1
                       
                PercT = TranscriptsOnTheDiagonal / float(TotalTranscripts)
                yaxis = ZeroInCase1 / float(TotalTranscripts)
                xaxis = ZeroInCase2 / float(TotalTranscripts)
                if diff == 0:
                    IntT = Int_Diagonal / float(TotalTranscripts)
                if v == '_rv':
                    RVDia += PercT
                    RVX += xaxis
                    RVY += yaxis
                    rvs.append(PercT)
                else:
                    AVGDiagonal += PercT
                    AVGYaxis += yaxis
                    AVGXaxis += xaxis
                    diagonals.append(PercT)
                    if v == "_s1":
                        AVGDiagonal_s1 += PercT
                    if v == "_s2":
                        AVGDiagonal_s2 += PercT
                    if v == "_s3":
                        AVGDiagonal_s3 += PercT
                   
                if diff == 0:
                    fi.write(ID + "\t" + v + "\t" + str(PercT) + "\t" + str(xaxis) + "\t" + str(yaxis) + "\t" + str(IntT)+ "\n")
                else:
                    fi.write(ID + "\t" + v + "\t" + str(PercT) + "\t" + str(xaxis) + "\t" + str(yaxis) + "\n")
               
                #Now we can do some comparisons
                plt.scatter(SharedOGTPMs, SharedModifiedTPMs, s=3)
                ax = plt.gca()
                #ax.axes.xaxis.set_ticklabels([])
                #ax.axes.yaxis.set_ticklabels([])
                plt.xlabel("Original Version log TPM", size = 8)
                if v == '_rv':
                    version = 'Reversed Complement'
                elif v == '_s1':
                    version = 'Shuffled Version 1'
                elif v == '_s2' :
                    version = 'Shuffled Version 2'
                else:
                    version = 'Shuffled Version 3'
                plt.title(key+" TPM Correlation between Original Version and " + version, size = 8, wrap=True)
                plt.ylabel(version+" log TPM", size = 8)
                plt.savefig(output_dir+key+'_'+ID + v+'_tpm_log.png',dpi=199)
                MDLines.clear()
                SharedOGTPMs.clear()
                SharedModifiedTPMs.clear()
                plt.clf()

       
        AVGDiagonal = float(AVGDiagonal) / 9
        AVGXaxis = float(AVGXaxis) / 9
        AVGYaxis = float(AVGYaxis) / 9
       
        AVGDiagonal_s1 = float(AVGDiagonal_s1) / 3
        AVGDiagonal_s2 = float(AVGDiagonal_s2) / 3
        AVGDiagonal_s3 = float(AVGDiagonal_s3) / 3
        fi.write("--------------------------------------------------------------------------\n")
        fi.write("Shuffled Average:\t" + str(AVGDiagonal) + "\t" + str(AVGXaxis) + "\t" + str(AVGYaxis)+"\n")
        fi.write("Shuffled 1 Average:\t" + str(AVGDiagonal_s1) + "\n")
        fi.write("Shuffled 2 Average:\t" + str(AVGDiagonal_s2) + "\n")
        fi.write("Shuffled 3 Average:\t" + str(AVGDiagonal_s3) + "\n")
        dia_np = np.array(diagonals)
        STD_shuf = np.std(dia_np)
        fi.write("Shuffled STD:\t" + str(STD_shuf)+"\n")
        RVDia = float(RVDia) / 3
        RVX = float(RVX) / 3
        RVY = float(RVY) / 3
        fi.write("RV Average:\t" + str(RVDia) + "\t" + str(RVX) + "\t" + str(RVY)+"\n")
        rv_np = np.array(rvs)
        STD_rv = np.std(rv_np)
        fi.write("RV STD:\t" + str(STD_rv)+"\n")
