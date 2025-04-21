
#create dictionary
#detect the codons, restore into dictionary
#if detect the terminals, stop
#AUG UUG GUG起始
#UAA UAG UGA终止

bases='ACGU'
codon_list=[a+b+c for a in bases for b in bases for c in bases]

gene_sequence='AUGUAUUAUUAUUAUACUACUUCAAUCUCGUCGUAGUAGCGUCGUUAGUAGUAUUAUUAUUAUUAAACUACUACUACUACUACU'
terminal=['UAA','UAG','UGA']
codon_dict={}
for codons in codon_list:
    codon_dict.setdefault(codons,0) 

i=0
while i>=0:
    codon=gene_sequence[3*i:3*i+3] 
    codon_dict[codon]+=1 
    i+=1
    if codon in terminal:
        break

max_codon = max(codon_dict, key=codon_dict.get)
print('The most frequent codon is',max_codon)

aa_dict={}
aa_dict['Alanine']=codon_dict['GCG']+codon_dict['GCA']+codon_dict['GCC']+codon_dict['GCU']
aa_dict['Cysteine']=codon_dict['UGC']+codon_dict['UGU']
aa_dict['Aspartic']=codon_dict['GAC']+codon_dict['GAU']
aa_dict['Glutamic_acid'] = codon_dict['GAG'] + codon_dict['GAA']
aa_dict['Phenylalanine'] = codon_dict['UUC'] + codon_dict['UUU']
aa_dict['Glycine'] = codon_dict['GGG'] + codon_dict['GGA'] + codon_dict['GGU'] + codon_dict['GGC']
aa_dict['Histidine'] = codon_dict['CAC'] + codon_dict['CAU']
aa_dict['Isoleucine'] = codon_dict['AUA'] + codon_dict['AUC'] + codon_dict['AUU']
aa_dict['Lysine'] = codon_dict['AAG'] + codon_dict['AAA']
aa_dict['Leucine'] = codon_dict['UUG'] + codon_dict['UUA'] + codon_dict['CUG'] + codon_dict['CUA'] + codon_dict['CUC'] + codon_dict['CUU']
aa_dict['Methionine'] = codon_dict['AUG']
aa_dict['Asparagine'] = codon_dict['AAC'] + codon_dict['AAU']
aa_dict['Proline'] = codon_dict['CCG'] + codon_dict['CCA'] + codon_dict['CCC'] + codon_dict['CCU']
aa_dict['Glutamine'] = codon_dict['CAG'] + codon_dict['CAA']
aa_dict['Arginine'] = codon_dict['CGG'] + codon_dict['CGA'] + codon_dict['CGC'] + codon_dict['CGU'] + codon_dict['AGG'] + codon_dict['AGA']
aa_dict['Serine'] = codon_dict['UCG'] + codon_dict['UCA'] + codon_dict['UCC'] + codon_dict['UCU'] + codon_dict['AGC'] + codon_dict['AGU']
aa_dict['Threonine'] = codon_dict['ACG'] + codon_dict['ACA'] + codon_dict['ACC'] + codon_dict['ACU']
aa_dict['Valine'] = codon_dict['GUA'] + codon_dict['GUC'] + codon_dict['GUU'] + codon_dict['GUG']
aa_dict['Tryptophan'] = codon_dict['UGG']
aa_dict['Tyrosine'] = codon_dict['UAC'] + codon_dict['UAU']

max_amino_acid = max(aa_dict, key=aa_dict.get)
print('The most frequent amino acid is',max_amino_acid)