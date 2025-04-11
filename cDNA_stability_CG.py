# check mRNA
# change mRNA -> cDNA:
#      find U and change into T
# measure CG propotion : moving average ,window size:100?
#    store information (series1)
# visualize CG propotion trend use heatmap


import re

# read or input
mRNA='ACGU'
# check mRNA
def mRNA_check(mRNA):
    pass

# change into cDNA
def cDNA(mRNA):
    mRNA=re.sub('U','T',mRNA)
    return mRNA

ws=100

# moving average
temL=[]
Results=[]
for i in range(0,len(mRNA)-ws):
    temL=mRNA[i:i+ws]
    temF=re.findall[r'[CG]',mRNA]
    temR=len(temF)
    Results.append(temR)
    temL,temR=[],0

