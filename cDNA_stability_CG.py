# check mRNA
# change mRNA -> cDNA:
#      find U and change into T
# measure CG propotion : moving average ,window size:100?
#    store information (series1)
# visualize CG propotion trend use heatmap

import matplotlib.pyplot as plt
import re
import seaborn as sns
import numpy as np
# read or input
mRNA='ACGAUCGAUCGCCCCCUGAUCGAUGGGGGGGGGGGGAUCGUAGCUACGUAGUAGUGUGUGUAGCUAGCUAGCAGCGAGCAACACAACGUUUUGUACAGUCGUAGUGUAGCUAGCUGUGUACUAGUCGAUGUGAUCGAUGCAUCGAUCGACUG'
# check mRNA, if not exit and report error

RNA_check=re.findall('[^AUCG]',mRNA)
if RNA_check:
    print('ERROR:This is not a RNA sequence')
    exit()

# change into cDNA
def cDNA(mRNA):
    mRNA=re.sub('U','T',mRNA)
    return mRNA

ws=50

# moving average
temL=[]
Results=[]
for i in range(0,len(mRNA)-ws):
    temL=mRNA[i:i+ws]
    temF=re.findall(r'[CG]',temL)
    temR=len(temF)
    Results.append(temR)
    temL,temR='',0

meanvalue=sum(Results)/len(Results)
#show heat map
data=np.array([Results])
plt.figure(figsize=(15, 3))  
sns.heatmap(data, cmap='coolwarm', vmin=0, vmax=ws, cbar_kws={'label': 'CG Content (%)'})
plt.title('CG Proportion Along cDNA Sequence (Sliding Window='+str(ws)+'bp)')
plt.xlabel('Window Position')
plt.yticks([]) 
plt.show()


# show line map
plt.figure(figsize=(10, 6))
sns.lineplot(data=Results)
plt.axhline(y=meanvalue, color='red', linestyle='--', label='Mean')
plt.title("CG Proportion Along cDNA")
plt.xlabel("position")
plt.ylabel("CG proportion")
plt.show()
