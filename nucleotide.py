
#create dictionary
#detect the codons, restore into dictionary
#if detect the terminals, stop
#AUG
#UAA UAG UGA

gene_sequence='AUGUAUUAUUAUUAUACUACUUCAAUCUCGUCGUAGUAGCGUCGUUAGUAGUAUUAUUAUUAUUAAACUACUACUACUACUACU'
terminal=['UAA','UAG','UGA']
codon_dict={}
i=0
while i>=0:
    codon=gene_sequence[3*i:3*i+3]
    print(codon)
    codon_dict.setdefault(codon, 0)
    codon_dict[codon]+=1
    i+=1
    if codon in terminal:
        break
print(codon_dict)