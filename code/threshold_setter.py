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
tester = 100

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

def contains_zero(pwm):
    """
    Detects the presence of a zero in a given PWM. Returns True or False.
    """

    for i in list(pwm.values()):
        if 0 in i:
            return True
    return False
        
def calculate_pssm(chosen_motif, pseudocounts=0.0):
    """
    chosen_motif: Chosen motif sequence to be used in creating the pwm
    pseudocounts: Default set to 0.0, otherwise added to the counts before calculating position-weight matrix.
    
    Returns the pssm (PositionSpecificScoringMatrix)
    """
    pwm = chosen_motif.counts.normalize(pseudocounts)
    if contains_zero(pwm):
        pwm = chosen_motif.counts.normalize(pseudocounts=0.1)
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

def getNegative(pos_seq):
    """
    pos_seq = dna sequence in the positive direction reading from the file
    
    Returns the negative counterpart of the positive sequence.
    """
    dict = {"A":'T','T':'A','G':'C','C':'G','-':'-'}
    negative = ""
    last_index = len(pos_seq) - 1
    while last_index > -1:
        negative += dict[pos_seq[last_index].upper()]
        last_index -= 1
    return negative

def sequence(ungapped, position, length):
    """
    Given an ungapped sequence and a positive or negative number (position),
    return the nucleotide at that position plus [length] nucleotides in the
    positive direction.
    """
    if position >= 0:
        return str(Seq(str(ungapped.seq), IUPAC.IUPACUnambiguousDNA())[position:position+length])
    else:
        return getNegative(str(Seq(str(ungapped.seq), IUPAC.IUPACUnambiguousDNA())[position:position+length]))


def positions(raw_sequence, cast_sequences, motif, ungapped, chosen_precision=10**4):
    """
    raw_sequence: A list of the sequences that have been ungapped but not casted.
    cast_sequences: A list of sequences that have been ungapped and casted.
    motif: the chosen motif file read in using read_motif
    precision: default precision set to 10**4 to be used in the pssm distribution.
    
    Returns a list of positions
    """
    pssm = calculate_pssm(motif)
    
    position_list = []
    len_and_ids = extract_len_id(raw_sequence, cast_sequences)
    for i in range(len(len_and_ids)):
        for position, score in pssm.search(cast_sequences[i], threshold = -50):
            pos = {'species': len_and_ids[i].get('species'), 'score':score, 
             'position': position, 'seq_len': len_and_ids[i].get('seq_len'),
             'nucleotide': sequence(ungapped[i], position, motif.length)}
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

def stand_cutoff(motif_path):
    """
    motif_path: (string) path of the 1 motif that you want to use (e.g. '../data/pwm/bcd_FlyReg.fm')
    Returns the 95th percentile of the scores, and the mean of the scores 
    """
    ofile = open("../data/alignments/random_fasta.txt", "w")
    for i in range(40):
        ofile.write(">" + generate_one_seq() + "\n")
    ofile.close()

    file = glob('../data/alignments/*random_fasta.txt')[0]

    curr_file = read_records(file)
    curr_motif = read_motif(motif_path)
    curr_raw = raw_sequence_ungap(curr_file)
    curr_cast = cast_sequence(curr_raw)

    pssm = calculate_pssm(curr_motif)

    position_list = []
    for i in (np.random.randint(0, len(curr_raw), size=5)):
    	for ps in pssm.search(curr_cast[i], threshold=0):
    		position_list.append(ps[1])

    percentile_95 = np.percentile(position_list, 95)
    return percentile_95


### Data Extraction
all_motifs = glob('../data/pwm/*.fm')

def filter_95_percentile(TFBS_df, stand_cutoff):
    """
    TFBS_df: A pandas df with the raw/aligned position, score, and file names.
    stand_cutoff: (integer) the 95% percentile of a randomly generated sequence;
        used to set the minimum score you want to keep in TFBS_df table.
        Generate stand_cutoff using the stand_cutoff function provided.
    
    returns the TFBS_df with only values above the 95th percentile.
    """
    minimum = stand_cutoff
    filtered_df = TFBS_df[TFBS_df["score"] >= minimum]
    return filtered_df
