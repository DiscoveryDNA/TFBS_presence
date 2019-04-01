from threshold_setter import *

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
        raw_df = pd.DataFrame(positions(curr_raw, curr_cast, curr_motif, curr_file))
        temp_df = positive_positions(raw_df)
        temp_df = define_sign(temp_df)

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
    
    upper_fifth = filter_95_percentile(TFBS_df, "../data/pwm/cad_FlyReg.fm")
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
def prelim_pipeline(alignment_file, motif_path, stand_cutoff=-50):

	table = calculate_one_dfs_TFBS(alignment_file, [motif_path])
	filtered_table = filter_95_percentile(table, stand_cutoff).drop(['seq_len', 'position'], axis=1)
	if (filtered_table.shape[0] == 0):
		print("No TFBS's detected in this alignment file.")
		return None
	else:
		return filtered_table

def pipeline(alignment_file, motif_path, stand_cutoff=0, percentile=80):
    """
    alignment_file: a string specifying the alignment file you want to use.
    motif_path: a string specifying the path to pwm file you want to use.
    stand_cutoff: 
    percentile: cutoff for high frequency of combination
    
    Returns:
    1) Table with unique alignment position, combinations at positions, and frequencies
    2) Table of high-frequency combinations, with high frequency defined as >= ? percentile
    3) Table containing information about the most frequency combination
    """
    # Filtered Table
    table = calculate_one_dfs_TFBS(alignment_file, [motif_path])
    filtered_table = filter_95_percentile(table, stand_cutoff).drop(['seq_len', 'position'], axis=1)
    if (filtered_table.shape[0] == 0):
        print("No TFBS's detected in this alignment file.")
        return None

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
    if results == None:
        return "No TFBS's detected in this alignment file."
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