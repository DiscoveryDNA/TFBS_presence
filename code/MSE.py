import numpy as np
import pandas as pd
import re
import Bio

import string
import random
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio import pairwise2
from IPython.display import Image
from Bio.pairwise2 import format_alignment
from glob import glob
import seaborn as sns
import matplotlib.pyplot as plt
from IPython.display import display
#get_ipython().magic('matplotlib inline')
plt.style.use('fivethirtyeight')

import os, sys
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna, generic_protein


def read_records(file):
    """
    file: a single .fa file within the directory.
    
    Returns the records from a .fa file as a list .
    """
    return list(SeqIO.parse(file, "fasta"))

def read_motif(motif):
    """
    motif: a single .fm file within the directory.
    
    Returns a motif.jaspar.Motif (the motif matrix) from the .fm file.
    """
    return motifs.read(open(motif), "pfm")


# The following functions clean the read in alligned sequences, taking out the '-' in them and then casting them as needed.
# We then feed this cleaned format to create a pssm (PositionSpecificScoringMatrix) with a patser_threshold applied on it
# to output a score with a corresponding index. Depending on the sign of the position we labeled the DNA strand as positive
# or negative and then saved all the calculated values as columns in a pandas df.

def raw_sequence_ungap(records):
    """
    records: a list of species within a single .fa file.
    
    Returns a list of Sequences with an ungapped('-') sequence.
    """
    raw_sequences = []
    for record in records:
        raw_sequences.append(SeqRecord(record.seq.ungap("-"), id = record.id))
    return raw_sequences

def cast_sequence(ungapped_sequence):
    """
    ungapped_sequence: a list with the sequence and id for all the species in a file.
    
    Returns a list sequences with the type cast as c.
    """
    casted = []
    for record in ungapped_sequence:
        casted.append(Seq(str(record.seq), IUPAC.IUPACUnambiguousDNA()))
    return casted

def calculate_pssm(chosen_motif, pseudocounts=0.0):
    """
    chosen_motif: Chosen motif sequence to be used in creating the pwm
    pseudocounts: Default set to 0.0, otherwise added to the counts before calculating position-weight matrix.
    
    Returns the pssm (PositionSpecificScoringMatrix)
    """
    pwm = chosen_motif.counts.normalize(pseudocounts)
    return pwm.log_odds()


def extract_len_id(raw_seq, raw_seq2):
    """
    raw_seq: A list of raw sequences that hasn't yet been cast.
    raw_seq2: A list of cast & ungapped sequences.
    
    Returns a list of dictionary pairs for the species and seq_len of each sequence.
    """
    raw_id = []
    record_length = []
    for seq in raw_seq:
        try:
            raw_id.append(seq.id)
        except AttributeError:
            raw_id.append('n/a')
            
    for record in raw_seq2:
            record_length.append(len(record))
    return [{'species': x, 'seq_len': y} for x, y in zip(raw_id, record_length)]


# Below are functions used to find the corresponding positions in the sequence and then mark them as
# negative or positive strands based on if the value was below or above 0 respectively. 


def positions(raw_sequence, cast_sequences, motif, chosen_precision=10**4):
    """
    raw_sequence: A list of the sequences that have been ungapped but not casted.
    cast_sequences: A list of sequences that have been ungapped and casted.
    motif: the chosen motif file read in using read_motif
    precision: default precision set to 10**4 to be used in the pssm distribution.
    
    Returns a list of positions
    """
    pssm = calculate_pssm(motif)
    distribution = pssm.distribution(background=motif.background, precision= chosen_precision)
    patser_threshold = distribution.threshold_patser()
    
    position_list = []
    len_and_ids = extract_len_id(raw_sequence, cast_sequences)
    for i in range(len(len_and_ids)):
        for position, score in pssm.search(cast_sequences[i], threshold= patser_threshold):
            pos = {'species': len_and_ids[i].get('species'), 'score':score, 
                         'position':position, 'seq_len': len_and_ids[i].get('seq_len') }
            position_list.append(pos)
    return position_list

def positive_positions(df):
    """
    df: The df with position, score, seq_len, and species.
    
    Returns a df where the position arguments are translated into positive ints.
    """
    temp_pos = df[df["position"] >= 0].copy()
    temp_neg = df[df["position"] < 0].copy()
    temp_pos["raw_position"] = temp_pos["position"]
    temp_neg["raw_position"] = temp_neg["seq_len"] + temp_neg["position"]  
    temp_together = temp_pos.append(temp_neg).reset_index().sort_values("index")

    return temp_together.set_index("index")

def define_sign(df):
    """
    df: The df with position, score, seq_len, species, and raw_position.
    
    Returns the df with a strand sign column appended on.
    """
    df["strand"] = np.where(df['position'] >= 0, 'positive', 'negative')
    return df

def merge_align_df(raw_df, aligned_seq, raw_id_len):
    """
    raw_df: A pandas df that contains the positions, score, seq_len, species, and strand orientation.
    aligned_seq: Original inputted sequence with the '-' still included.
    raw_id_len: The length of the o
    
    Returns a df with alignned index appended to raw_df sorted by ['species', 'raw_position'] 
    with N/A values dropped.
    
    """
    remap_list = []
    nuc_list = ['A', 'a', 'G', 'g', 'C', 'c', 'T', 't', 'N', 'n']

    for i in range(raw_id_len):
        counter = 0
        for xInd, x in enumerate(aligned_seq[i].seq):    
            if x in nuc_list:
                remaps = {'raw_position': counter, 'align_position':xInd, 'species':aligned_seq[i].id}
                counter += 1
                remap_list.append(remaps)
            
    remap_DF = pd.DataFrame(remap_list)
    TFBS_map_DF_all = pd.merge(raw_df, remap_DF, on=['species', 'raw_position'], how='outer')
    TFBS_map_DF_all = TFBS_map_DF_all.sort_values(by=['species','align_position'], ascending=[True, True])

    return TFBS_map_DF_all.dropna()   

def save_df(TFBS_df, align_file, motif_file):
    """
    TFBS_df: A pandas df with the raw/aligned position, score, and file names.
    align_file: Name of read in alignment file.
    motif_file: Name of read in motif file.
    
    Returns a csv copy of the df stored in the current directory
    """
    return TFBS_df.to_csv(motif_file + ".csv", sep='\t', na_rep="NA")


# # Creating a Standard Score Distribution
# The following 4 cells generate a random sequence, finds the 95 percentile cutoff point in the scores of the random sequence (differs depending on which motif you use), and uses this as the standard cutoff point for your own observed scores.


# Randomly generates a sequence
def generate_one_seq():
    header = "VT0000|1|MEMB001A|-|"+str(np.random.randint(1000,9000))
    random_list = [random.choice("ACGT") for i in range(5000)]
    sequence = "".join(random_list) # 1 long string
    return header + "\n" + sequence

def generate_fasta(path):
    """
    path: a string indicating the location in which random file will be created / used
    """
    ofile = open(path, "w")
    for i in range(40):
        ofile.write(">" + generate_one_seq() + "\n")
    ofile.close()
    return glob(path)[0]

def stand_cutoff(motifs, motif_to_use):
    """
    files: one fasta file (string)
    motifs: a list of motifs
    motif_to_use: a number to index into (chose 1 motif you want to use)
    Returns the 95th percentile of the scores, and the mean of the scores 
    """
    ofile = open("../data/alignments/random_fasta.txt", "w")
    for i in range(40):
        ofile.write(">" + generate_one_seq() + "\n")
    ofile.close()

    file = glob('../data/alignments/*random_fasta.txt')[0]

    motif = motifs[motif_to_use]
    curr_file = read_records(file)
    curr_motif = read_motif(motif)
    curr_raw = raw_sequence_ungap(curr_file)
    curr_cast = cast_sequence(curr_raw)

    len_raw = len(extract_len_id(curr_raw, curr_cast))
    raw_df = pd.DataFrame(positions(curr_raw, curr_cast, curr_motif))
    temp_df = positive_positions(raw_df)
    percentile_95 = np.percentile(temp_df["score"], 95)
    standard_mean = temp_df["score"].mean()
    return (percentile_95, standard_mean)


### Data Extraction
all_motifs = glob('../data/pwm/*.fm')
cutoff_dict = {}
cutoff_dict["bcd_FlyReg.fm"] = stand_cutoff(all_motifs, 0)[0]
cutoff_dict["cad_FlyReg.fm"] = stand_cutoff(all_motifs, 1)[0]
cutoff_dict["eve_new6.fm"] = stand_cutoff(all_motifs, 2)[0]
cutoff_dict["zelda_.fm"] = stand_cutoff(all_motifs, 3)[0]

def filter_95_percentile(TFBS_df, motif_used):
    """
    motif_used: a string that specifies the motif you used (e.g. "cad_FlyReg.fm")
    TFBS_df: A pandas df with the raw/aligned position, score, and file names.
    
    returns the TFBS_df with only values above the 95th percentile.
    """
    minimum = cutoff_dict[motif_used]
    filtered_df = TFBS_df[TFBS_df["score"] >= minimum]
    return filtered_df

def calculate_one_dfs_TFBS(file, all_motifs):
    """
    file: A single .fa file from the current directory.
    all_motifs: A list of all the .fm files in the current directory.
    
    Returns a df with the positions, score, and orientation saved into the current directory for all files.
    """

    for motif in all_motifs:
        curr_file = read_records(file)
        curr_motif = read_motif(motif)
        curr_raw = raw_sequence_ungap(curr_file)
        curr_cast = cast_sequence(curr_raw)
        curr_motif_name = os.path.basename(motif) 

        len_raw = len(extract_len_id(curr_raw, curr_cast))
        raw_df = pd.DataFrame(positions(curr_raw, curr_cast, curr_motif))
        temp_df = positive_positions(raw_df)
        temp_df = define_sign(temp_df)

        align_name = re.split(r'_', file)[-1]
        final_df = merge_align_df(temp_df, curr_file, len_raw)
        final_df['motif'] = curr_motif_name

    return final_df

def map_raw_seq(TFBS_df, seq_file):
    """
    TFBS_df: The TFBS_df with entries only above the 95th percentile.
    seq_file: a single .fa file within the directory.
    
    returns a df with scores in the upper 5th percentile along with the species'
    corresponding raw sequence.
    """
    
    values = [[record.seq, record.description] for record in seq_file]
    sequences = pd.DataFrame(values)
    sequences.rename(columns={0 : 'Sequence', 1: 'Description'}, inplace=True)
    sequences = sequences.set_index("Description")
    mapping = sequences.to_dict("index")
    
    upper_fifth = filter_95_percentile(TFBS_df, "cad_FlyReg.fm")
    upper_fifth = upper_fifth[["species", "raw_position",'align_position', "strand"]]
    upper_fifth["sequence"] = upper_fifth["species"].apply(lambda x: str(mapping.get(x).get("Sequence")))
    
    upper_fifth = upper_fifth.reset_index().drop("index", axis=1)
    
    return upper_fifth

def positive_index_seq(matched_df, motif):
    """
    matched_df: df with the species matched with their corresponding raw sequence.
    motif: a single .fm file within the directory.
    
    returns the df with appended columns of the DNA sequences/indices starting from the 
    raw_index from positively oriented scores and then going on until the proper number
    of letters have been found.
    """
    motif = str(first_motif.consensus)
    motif_len = len(motif)
    
    extracted_seqs = []
    extracted_indices = []
    
    pos = matched_df[matched_df["strand"] == "positive"]
    final = pos.copy()
    for i in np.arange(pos.shape[0]):
        curr_row = pos.iloc[i]
        curr_seq = curr_row["sequence"]
        index = int(curr_row["raw_position"])
        
        selected_seq = curr_seq[index:index + motif_len]     
        letters = re.findall(r'[^-]', selected_seq)
        additional_indicies = 0
        
        if len(letters) != motif_len: 
            while len(re.findall(r'[^-]', selected_seq)) != motif_len:
                additional_indicies += 1
                selected_seq = curr_seq[index:index + motif_len + additional_indicies]
                 
        selected_index = "{0}:{1}".format(index, index + motif_len + additional_indicies)
        extracted_seqs.append(selected_seq)
        extracted_indices.append(selected_index)
        
    final["extracted_seq"] = extracted_seqs
    final["extracted_indices"] = extracted_indices
    return final

def negative_index_seq(matched_df, motif):
    """
    matched_df: df with the species matched with their corresponding raw sequence.
    motif: a single .fm file within the directory.
    
    returns the df with appended columns of the DNA sequences/indices starting from the 
    raw_index from negatively oriented scores and then going on until the proper number
    of letters have been found.
    """
    motif = str(first_motif.consensus)
    motif_len = len(motif)
    
    extracted_seqs = []
    extracted_indices = []
    
    d={'A':'T','T':'A', 'C':'G', 'G':'C'}
    
    neg = matched_df[matched_df["strand"] == "negative"]
    final = neg.copy()
    for i in np.arange(neg.shape[0]):
        curr_row = neg.iloc[i]
        curr_seq = curr_row["sequence"]
        
        index = int(curr_row["raw_position"])
        
        selected_seq = curr_seq[-(index + motif_len + 1): -(index + 1)]     
        letters = re.findall(r'[^-]', selected_seq)
        additional_indicies = 0
        
        if len(letters) != motif_len: 
            while len(re.findall(r'[^-]', selected_seq)) != motif_len:
                additional_indicies += 1
                selected_seq = curr_seq[-(index + motif_len + additional_indicies + 1): -(index + 1)]
                 
        selected_index = "{0}:{1}".format(-(index + motif_len + additional_indicies + 1), -(index + 1)) 
        
        swapped_seq = ''.join(d[s] if s in d else s for s in selected_seq)
        reversed_seq = swapped_seq[::-1]
        extracted_seqs.append(reversed_seq)
        extracted_indices.append(selected_index)
        
    final["extracted_seq"] = extracted_seqs
    final["extracted_indices"] = extracted_indices
    return final


def find_motif_pipeline(TFBS_df, seq_file, motif):
    """
    TFBS_df: A pandas df with the raw/aligned position, score, and file names.
    seq_file: a single .fa file within the directory.
    motif: a single .fm file within the directory.
    
    returns the execution of the pipeline for filtering out the calculated motif locations
    and the actual DNA sequence.
    """
    matched = map_raw_seq(TFBS_df, seq_file)
    final_p = positive_index_seq(matched, motif)
    final_n = negative_index_seq(matched, motif)
    final = final_p.reset_index().merge(final_n.reset_index(), how='outer').sort_values('index').set_index('index')
    final["gap_less_seq"] = [re.sub(r'-', '', x) for x in final["extracted_seq"]]
    return final[["species", "raw_position", "strand","gap_less_seq", "extracted_indices", 'align_position']]


def calculate_one_TFBS(file, motif):
    """
    file: A single .fa file in the current directory.
    motif: A single .fm file in the current directory.
    
    Returns a csv with the positions, score, and orientation saved into the current directory for the specified 
    .fa and .fm files.
    """
    curr_file = read_records(file)
    curr_motif = read_motif(motif)
    curr_raw = raw_sequence_ungap(curr_file)
    curr_cast = cast_sequence(curr_raw)
            
    len_raw = len(extract_len_id(curr_raw, curr_cast))    
    raw_df = pd.DataFrame(positions(curr_raw, curr_cast, curr_motif))
    temp_df = positive_positions(raw_df)
    temp_df = define_sign(temp_df)
            
    align_name = re.split(r'_', files[0])[-1]
    final_df = merge_align_df(temp_df, curr_file, len_raw)
    save_df(final_df, align_name, motif)
    
    new_path = curr_motif.weblogo("{}.png".format(motif), format= 'PNG', stack_width= 'large')
    display(Image("{}.png".format(motif), width= 200, height= 300))

    indexed_df = find_motif_pipeline(final_df, curr_file, curr_motif)
    testing_pos = indexed_df[indexed_df["strand"] == 'positive']
    testing_neg = indexed_df[indexed_df['strand'] == 'negative']

    fig = plt.figure(figsize=(20, 15))
    plt.subplot(2, 1, 1)
    sns.stripplot(x="align_position", y="species", hue="strand", data=testing_pos, dodge=True, linewidth=1)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,title='Orientation', fontsize='x-large')
    plt.xlabel('')
    plt.title('{} & {} Alignment'.format(align_name, motif))

    plt.subplot(2, 1, 2)
    sns.stripplot(x="align_position", y="species", hue="strand", data=testing_neg,
                  palette="Set2", dodge=True, linewidth=1)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,title='Orientation', fontsize='x-large')
    plt.xlabel('')
    
    #fig.savefig('{} & {}.png'.format(align_name, motif), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.5)
    return final_df


def calculate_all_TFBS(files, all_motifs):
    """
    files: A list of all the .fa files in the current directory.
    all_motifs: A list of all the .fm files in the current directory.
    
    Returns a csv with the positions, score, and orientation saved into the current directory for all files.
    """
    for file in files:
        for motif in all_motifs:
            calculate_one_TFBS(file, motif)

def frequent_patterns(table, percentile=80):
    """
    alignment_num: integer specifying which alignment to analyze
    tables: a list of tables like group_by_combo_tables
    percentile: (integer) limit frequencies to only above this percentile
    
    Returns:
    A table, like those in group_by_combo_tables, except w/only the largest counts
    """
    counts = list(table["Count"])
    return table[table["Count"] >= np.percentile(counts, percentile)]

def most_frequent_pattern(table):
    """
    alignment_num: integer specifying which alignment to analyze
    tables: a list of tables like group_by_combo_tables
    
    Returns:
    A single-row table with the TFBS combo with the LARGEST count in the sequence
    """
    max_count = max(list(table["Count"]))
    return table[table["Count"] == max_count]

def frequent_across_alignments(align_list, motif_path, percentile):
    """
    align_list: a list of multiple alignment files
    motif_path: a string - the path to pwm file you want to use.
    percentile: cutoff for high frequency of combination
    """
    freq_tables = []
    for i in align_list:
        fp = pipeline(i, motif_path, percentile)[2]
        short_name = i[48:]
        fp["Alignment File"] = [short_name] * fp.count()[0]
        freq_tables.append(fp)
    
    frequent_all_alignments = pd.concat(freq_tables).sort_values(by=['Count'], ascending=False)
    return frequent_all_alignments

### Wrapping Into One Package

def pipeline(alignment_file, motif_path, percentile=80):
    """
    alignment_file: a string specifying the alignment file you want to use.
    motif_path: a string specifying the path to pwm file you want to use.
    percentile: cutoff for high frequency of combination
    
    Returns:
    1) Table with unique alignment position, combinations at positions, and frequencies
    2) Table of high-frequency combinations, with high frequency defined as >= ? percentile
    3) Table containing information about the most frequency combination
    """
    # Filtered Table
    table = calculate_one_dfs_TFBS(alignment_file, [motif_path])
    filtered_table = filter_95_percentile(table, "cad_FlyReg.fm").drop(['seq_len', 'position'], axis=1)

    # Position_species_table
    positions_with_TFBS = list(filtered_table["align_position"])
    combos = {} # Dictionary with position index : list of species that have TFBS at this position
    for i in positions_with_TFBS:
        combos[i] = list(filtered_table[filtered_table["align_position"] == i]["species"])
    d = {'Alignment Positions': list(combos.keys()),
         'Species with TFBS at Position': list(combos.values())}
    positions_species_table = pd.DataFrame(data=d)
    positions_species_table["Species with TFBS at Position"] = [tuple(i) for i in positions_species_table["Species with TFBS at Position"]]

    # Group_by_combo
    group_by_combo = positions_species_table.groupby("Species with TFBS at Position",sort=False).count().reset_index()
    group_by_combo.columns = ["Species with TFBS at Position", "Count"]
    
    freq_patterns = frequent_patterns(group_by_combo, percentile)
    most_freq_pattern = most_frequent_pattern(group_by_combo)

    return filtered_table, group_by_combo, freq_patterns, most_freq_pattern

def view_data(results, view = "all"):
    """
    results: the list of tables outputted from the pipline function.

    view: a string ("filtered", "grouped", "frequent", or "most frequent") specifying the type of result you want to view.
        "filtered": a table containing the score, species, raw position, strand, align_position, and motif of each TFBS.
        "grouped": a table containing the 1) Species with TFBS at Position, and 2) Count of that combination
        "frequent": a table, like grouped, except that it includes only those with large counts
        "most frequent": a single row table with the most frequently occurring combination and its count.
    
    Provides an easier and clearer way to view specific parts of the output you want to see.
    """
    if view == "filtered":
        return results[0]
    elif view == "grouped":
        return results[1]
    elif view == "frequent":
        return results[2]
    elif view == "most frequent":
        return results[3]
    else:
        print("request not recognized")


# Deleted from original MSE:
# * pwm_threshold
# * calculate_dfs_TFBS
