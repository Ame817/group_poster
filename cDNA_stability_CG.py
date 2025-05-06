# check mRNA
# change mRNA -> cDNA:
#      find U and change into T
# 设置窗口大小：50
# measure CG propotion : 
#     从第一个开始：计算【1到50】50个base中CG的百分比
#     将窗口向右滑动一个base，计算从【2到51】50个base中CG的百分比
#     以此类推，一直往右边推，推到最右边为止
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
def cDNA_turn(mRNA):
    cDNA=re.sub('U','T',mRNA)
    return cDNA

cDNA=cDNA_turn(mRNA)

ws=50

# moving average
temL=[]     #临时存放50个base的列表
Results=[]  #存放结果的list
for i in range(0,len(cDNA)-ws): # 从第一个开始，一直往右边推
    temL=cDNA[i:i+ws]  
    temF=re.findall(r'[CG]',temL)
    temR=len(temF)
    Results.append(temR) #记录多少个CG
    temL,temR='',0
Results=[t/ws for t in Results] #把CG换算成比例
#show heat map
data=np.array([Results]) 
plt.figure(figsize=(15, 3))  
sns.heatmap(data, cmap='coolwarm', vmin=0, vmax=1, cbar_kws={'label': 'CG Content (%)'})
plt.title('CG Proportion Along cDNA Sequence (Sliding Window='+str(ws)+'bp)')
plt.xlabel('Window Position')
plt.yticks([]) 
plt.show()

