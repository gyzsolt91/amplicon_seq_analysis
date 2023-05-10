import pandas as pd
from Bio import SeqIO
import os
from configs import RunConfig, FileConfig
from functions import make_frequency_table, save_seqlogo
import regex


def analysis(runconfig: RunConfig, fileconfig: FileConfig):
    samples = runconfig.samples
    barcodes = runconfig.barcodes.dict()
    possible_bases = runconfig.possible_bases
    in_dir = fileconfig.out_dir
    out_dir = fileconfig.out_dir
    figure_dir = fileconfig.figures_dir
    ref_seq = runconfig.ref_seq
    ref_seq_prefix = ref_seq[10:21]
    ref_seq_postfix = ref_seq[-16:-5]
    prefix_regex = f'({ref_seq_prefix})' + '{e<=1}'
    postfix_regex = f'({ref_seq_postfix})' + '{e<=1}'
    start = regex.search(prefix_regex, ref_seq).end()
    end = regex.search(postfix_regex, ref_seq).start()
    expected = len(ref_seq[start:end])

    os.makedirs(figure_dir, exist_ok=True)

    for sample in samples:
        print(f'Analyzing {sample} sample.')
        table_path = os.path.join(in_dir, f'{sample}_freq_table.csv')
        fig_path = os.path.join(figure_dir, f'{sample}_seqlogo.jpeg')
        input_file = os.path.join(in_dir, f'{sample}_matching_length.fastq')
        sequences = [rec.seq for rec in SeqIO.parse(input_file, "fastq")]
        unique_seq_path = os.path.join(out_dir, f'{sample}_unique_seqs.csv')
        mismatch_path = os.path.join(out_dir, f'{sample}_mismatch_combinations.csv')
        unique_sequences = {}
        mismatch_combination_dict = {
            'A': [0, 0, 0, 0, 0],
            'C': [0, 0, 0, 0, 0],
            'G': [0, 0, 0, 0, 0],
            'T': [0, 0, 0, 0, 0],
            'N': [0, 0, 0, 0, 0]
        }
        mismatch_combination_freq = pd.DataFrame(mismatch_combination_dict, index=mismatch_combination_dict.keys())

        for rec in SeqIO.parse(input_file, "fastq"):
            if rec.seq not in unique_sequences.keys():
                unique_sequences[rec.seq] = 1
            else:
                unique_sequences[rec.seq] += 1

            upstream_mismatch = rec.seq[9]
            downstream_mismatch = rec.seq[31]
            mismatch_combination_freq.loc[upstream_mismatch, downstream_mismatch] += 1

        unique_sequences_list = []
        for key, value in unique_sequences.items():
            unique_sequences_list.append({'sequence': key, 'count': value})

        unique_seq_df = pd.DataFrame(unique_sequences_list)
        unique_seq_df.to_csv(unique_seq_path, index=False)
        mismatch_combination_freq.to_csv(mismatch_path)

        frequency_df = make_frequency_table(sequences=sequences,
                                            length_of_seq=expected,
                                            possible_bases=possible_bases,
                                            path=table_path)

        #print(frequency_df)

        save_seqlogo(frequency_df=frequency_df,
                     sample=sample,
                     path=fig_path)

