import pandas as pd
from Bio import SeqIO, AlignIO, Seq
from collections import Counter
import numpy as np

def parse_uclust(infile,fasta,outfasta,outfile):
    header = ['Type','Cluster','Size','%Id','Strand','Qlo','Tlo','Alignment','Query','Target']
    uclust = pd.read_csv(infile, sep='\t', names=header, index_col=False)
    centroids = uclust[uclust['Type']=='C']
    c_list = centroids.iloc[:,8].values
    c_list = [c.split(" ")[0] for c in c_list]

    sequences = []
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in c_list:
            sequences.append(record)

    seq2 = []
    seq3 = []
    for i,record in enumerate(sequences):
        if not record.id in seq2: 
            seq2.append(record.id)
            seq3.append(record)
    with open(outfasta, "w") as output_handle:
        SeqIO.write(seq3, output_handle, "fasta")
    
    uclust.to_csv(outfile,index=False)

def clean_aln(infile,outfile):
    aln = AlignIO.read(infile,format='fasta')
    mask = []
    for i in range(0,aln.get_alignment_length()):
        mask.append(Counter(aln[:,i])['-']/len(aln[:,i])<0.95)

    np_aln = np.array(aln)
    np_aln = np_aln[:,mask]

    aln_faa = SeqIO.parse(infile,format='fasta')
    sequences = []
    for record,i in zip(aln_faa,range(0,len(aln))):
        record.seq = Seq.Seq("".join(np_aln[i,:]), Seq.Alphabet.SingleLetterAlphabet())
        sequences.append(record)

    with open(outfile, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
        
def add_type(type_file,seq_file,outfile):
    labels = pd.read_csv(type_file)
    uclust_data = pd.read_csv(seq_file)
    uclust_data = uclust_data[uclust_data['Type'] !='S']
    uclust_data.loc[uclust_data['Target'] == '*','Target'] = uclust_data.loc[uclust_data['Target'] == '*','Query']
    uclust_data = uclust_data.merge(labels, left_on='Query', right_on='ID',how='left')
    labeled_leaves = uclust_data[~pd.isna(uclust_data['type'])]
    color_map = {'I': '#28B463',
                 'II': '#E74C3C',
                 'II/III':'#AF7AC5',
                 'IIIa':'#AED6F1',
                 'IIIb':'#3498DB',
                 'IIIc':'#1F618D',
                 'IIIlike':'#5D6D7E',
                 'IV':'#F4D03F',
                 'IVlike':'#F8C471',
                 'unknown':'#F442D4'}
    lines = labeled_leaves['Target'].apply(lambda x: x.split(' ')[0]).values +[',label,node,'+color_map[x]+',1,normal\n' for x in labeled_leaves['type'].values]
    lines = set(lines)
    with open('../data/itol_legend_template.txt','r') as file:
        with open(outfile, "w") as f1:
            for row in file:
                f1.write(row)
            for line in lines:
                f1.write(line)
            file.close()
            f1.close()

def add_kinetic(kinetic_file,synth_file,seq_file,outfile):
    kinetic_data = pd.DataFrame([x.description for x in SeqIO.parse(kinetic_file, "fasta")],columns=['kinetic_ID'])
    uclust_data = pd.read_csv(seq_file)
    uclust_data = uclust_data[uclust_data['Type'] !='S']
    uclust_data.loc[uclust_data['Target'] == '*','Target'] = uclust_data.loc[uclust_data['Target'] == '*','Query']
    uclust_data = uclust_data.merge(kinetic_data, left_on='Query', right_on='kinetic_ID',how='left')
    
    synth_data = pd.DataFrame([x.description for x in SeqIO.parse(synth_file, "fasta")],columns=['syn_ID'])
    uclust_data = uclust_data.merge(synth_data, left_on='Query', right_on='syn_ID',how='left')
    
    uclust_data['kinetic_flag'] = '-1'
    uclust_data['syn_flag'] = '-1'

    kinetic_centroid = uclust_data.loc[~pd.isna(uclust_data['kinetic_ID']),'Target'].unique()
    syn_centroid = uclust_data.loc[~pd.isna(uclust_data['syn_ID']),'Target'].unique()
    
    uclust_data.loc[uclust_data['Target'].isin(kinetic_centroid),'kinetic_flag'] = '1'
    uclust_data.loc[uclust_data['Target'].isin(syn_centroid),'syn_flag'] = '1'

    lines = uclust_data['Target'].apply(lambda x: x.split(' ')[0]).values + ','+ uclust_data['kinetic_flag'].values+','+uclust_data['syn_flag'].values+'\n'
    unique_lines = np.unique(lines)
    
    with open('../data/kinetic_sampling_legend.txt','r') as file:
        with open(outfile, "w") as f1:
            for row in file:
                f1.write(row)
            for line in unique_lines:
                f1.write(line)
            file.close()
            f1.close()