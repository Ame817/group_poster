
#create dictionary
#detect the codons, restore into dictionary
#if detect the terminals, stop
#AUG UUG GUG initiator(assume that the RNA sequence start from AUG)
#UAA UAG UGA terminal


import matplotlib.pyplot as plt
import re

bases='ACGU'
codon_list=[a+b+c for a in bases for b in bases for c in bases]

#the input should be RNA sequence.这里放了DNA是因为我还没找到能直接用的RNA序列 后面要改
gene_sequence='ATGTTGAATAGTTCAAGAAAATATGCTTGTCGTTCCCTATTCAGACAAGCGAACGTCTCAATAAAAGGACTCTTTTATAATGGAGGCGCATATCGAAGAGGGTTTTCAACGGGATGTTGTTTGAGGAGTGATAACAAGGAAAGCCCAAGTGCAAGACAACCACTAGATAGGCTACAACTAGGTGATGAAATCAATGAACCAGAGCCTATTAGAACCAGGTTTTTTCAATTTTCCAGATGGAAGGCCACCATTGCTCTATTGTTGCTAAGTGGTGGGACGTATGCCTATTTATCAAGAAAAAGACGCTTGCTAGAAACTGAAAAGGAAGCAGATGCTAACAGAGCTTACGGTTCAGTAGCACTTGGCGGTCCTTTCAATTTAACAGATTTTAATGGTAAGCCTTTCACTGAGGAGAATTTGAAGGGTAAGTTTTCCATTTTATACTTTGGATTCAGTCATTGCCCCGACATTTGTCCAGAAGAGCTTGACAGATTAACGTATTGGATTTCTGAATTAGATGATAAAGACCATATAAAGATACAGCCATTGTTTATCTCATGTGATCCTGCAAGAGATACACCGGATGTCTTGAAAGAGTACTTAAGCGATTTTCACCCAGCTATCATTGGTTTAACCGGTACGTACGACCAAGTGAAAAGCGTATGCAAAAAATACAAGGTATATTTTTCAACTCCACGTGATGTCAAGCCCAACCAGGATTACTTAGTGGACCATTCGATATTTTTCTATTTGATCGACCCTGAAGGACAGTTTATCGATGCGTTGGGAAGAAACTACGATGAGCAATCTGGTCTCGAAAAGATTCGTGAACAAATTCAGGCGTATGTGCCAAAGGAAGAACGGGAGCGTAGGTCAAAAAAATGGTACTCTTTTATCTTCAATTGA'
RNA_sequence=gene_sequence.replace('T', 'U')

RNA_check=re.findall('[^AUCG]',RNA_sequence)
if RNA_check:
    print('ERROR:This is not a RNA sequence')
else:
    def frequent_codon(RNA_sequence,silent=False):    
        terminal=['UAA','UAG','UGA']
        codon_dict={}
        for codons in codon_list: # create condons index in dict and set defalt value
            codon_dict.setdefault(codons,0) 
        i=0
        while i>=0: # count condons
            codon=RNA_sequence[3*i:3*i+3] 
            codon_dict[codon]+=1 
            i+=1
            if codon in terminal: #if meet terminal condon, exit 
                break
        max_codon=max(codon_dict, key=codon_dict.get)
        if not silent: # if not scilence, output
            print('The most frequent codon in this sequence is',max_codon)
        return(codon_dict)

    def frequent_amino_acids(RNA_sequence,silent=False):
        codon_dict=frequent_codon(RNA_sequence,silent=True) # get condon dict
        aa_dict={}
        #count frequency
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

        max_amino_acid=max(aa_dict, key=aa_dict.get) #get most frequent amino acid
        if not silent:
            print('The most frequent amino acid is',max_amino_acid)
        return(aa_dict)

    def aminoacids_piechart(RNA_sequence): #generate pie chart, in decsending order.
        aa_dict = frequent_amino_acids(RNA_sequence, silent=True)
        labels = list(aa_dict.keys())
        values = list(aa_dict.values())
        sort_list = sorted(zip(labels, values), key=lambda x: x[1], reverse=False) # sort value in descending order
        sorted_labels, sorted_values = zip(*sort_list)
        cmap = plt.get_cmap('coolwarm') # set color in red to blue
        colors = [cmap(j / len(sorted_values)) for j in range(len(sorted_values))]
        explode = [0.05] * len(sorted_values)
    
        plt.figure(figsize=(8, 8))  
        plt.pie(sorted_values, labels=sorted_labels,autopct='%1.1f%%', startangle=140, colors=colors,explode=explode)
        plt.axis('equal')
        plt.title('The frequency of amino acids')
        plt.show()

frequent_codon(RNA_sequence)
frequent_amino_acids(RNA_sequence)
aminoacids_piechart(RNA_sequence)