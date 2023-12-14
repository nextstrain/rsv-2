import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

def sequence_to_int_array(s, fill_value=110, fill_gaps=True):
    seq = np.frombuffer(str(s).lower().encode('utf-8'), dtype=np.int8).copy()
    if fill_gaps:
        seq[(seq!=97) & (seq!=99) & (seq!=103) & (seq!=116)] = fill_value
    else:
        seq[(seq!=97) & (seq!=99) & (seq!=103) & (seq!=116) & (seq!=45)] = fill_value
    return seq

def get_similarity(alignment, reference):
    ref_array = sequence_to_int_array(SeqIO.read(reference,"fasta").seq, fill_gaps=False)
    similarity = {}
    for seq in SeqIO.parse(alignment,"fasta"):
        seq_array = sequence_to_int_array(seq.seq, fill_gaps=False)
        good_indices = seq_array!=45
        similarity[seq.id] = np.mean(seq_array[good_indices]==ref_array[good_indices])

    return similarity

if __name__=="__main__":
    dfs, a_list, b_list, a_sequences, b_sequences, listofdictionaries, a_metadata, b_metadata = ([] for i in range(8))
    to_a = ['a_1', 'a_2', 'a_3']
    to_b = ['b_1', 'b_2', 'b_3']
    seq = ["data/a/1_sequences.aligned.fasta", "data/a/2_sequences.aligned.fasta", "data/a/3_sequences.aligned.fasta",
           "data/b/1_sequences.aligned.fasta", "data/b/2_sequences.aligned.fasta", "data/b/3_sequences.aligned.fasta"]
    references = ["config/a_1_reference.fasta", "config/a_2_reference.fasta", "config/a_3_reference.fasta",
                  "config/b_1_reference.fasta", "config/b_2_reference.fasta", "config/b_3_reference.fasta"]
    accessions = [pd.read_csv("data/a/1_metadata.tsv", sep='\t').accession,
                  pd.read_csv("data/a/2_metadata.tsv", sep='\t').accession,
                  pd.read_csv("data/a/3_metadata.tsv", sep='\t').accession,
                  pd.read_csv("data/b/1_metadata.tsv", sep='\t').accession,
                  pd.read_csv("data/b/2_metadata.tsv", sep='\t').accession,
                  pd.read_csv("data/b/3_metadata.tsv", sep='\t').accession]
    everything = set().union(*accessions)

    all_similarities = defaultdict(list)
    for filename, reference in zip(seq, references):
        similarity = get_similarity(filename, reference)

        for accession in everything:
            all_similarities[accession].append(similarity.get(accession, 0.001))

    all_similarities_df = pd.DataFrame(all_similarities, index = to_a + to_b).T

    max_similarity = all_similarities_df.max(axis=1)

    all_similarities_df = all_similarities_df.loc[max_similarity>0.70,:]

    a_or_b = all_similarities_df.idxmax(axis=1).apply(lambda x: 'A' if x[0]=='a' else 'B')

    for record in SeqIO.parse("data/sequences.fasta", "fasta"):
        if record.id in a_or_b:
            if a_or_b[record.id]=='A':
                a_sequences.append(record)
            elif a_or_b[record.id]=='B':
                b_sequences.append(record)

    SeqIO.write(a_sequences, "data/a/sequences.fasta","fasta")
    SeqIO.write(b_sequences, "data/b/sequences.fasta", "fasta")

    metadata = pd.read_csv("data/metadata.tsv", sep="\t", index_col='accession')
    original_columns = metadata.columns
    metadata.drop_duplicates(keep='first', inplace=True)
    metadata['type'] = a_or_b

    a_metadata = pd.DataFrame(data=metadata.loc[metadata['type']=='A'],
                              columns=original_columns)
    b_metadata = pd.DataFrame(data=metadata.loc[metadata['type']=='B'],
                              columns=original_columns)
    a_metadata.to_csv('data/a/metadata_sorted.tsv', sep="\t")
    b_metadata.to_csv('data/b/metadata_sorted.tsv', sep="\t")
