import sys
from Bio import SeqIO
from Bio.Seq import Seq
import os
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from configs import RunConfig, FileConfig
from functions import make_frequency_table
import logomaker
import regex


def analysis(runconfig: RunConfig, fileconfig: FileConfig):
    barcodes = runconfig.barcodes.dict()
    possible_bases = runconfig.possible_bases
    input_dir = fileconfig.out_dir
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

    for bc in barcodes.keys():
        print(f'Analyzing {barcodes[bc]} sample.')
        filename = os.path.join(input_dir, f'{barcodes[bc]}_matching_length.fastq')
        sequences = [rec.seq for rec in SeqIO.parse(filename, "fastq")]

        frequency_df = make_frequency_table(sequences=sequences,
                                            length_of_seq=expected,
                                            possible_bases=possible_bases)

        print(frequency_df)
        print('Starts making seqlogo')
        logo = logomaker.Logo(frequency_df)
        plt.show()




    #
    #
    #
    #
    # barcodes = runconfig.barcodes.dict()
    # odir = fileconfig.out_dir
    # idir = fileconfig.out_dir
    # fdir = "figures"
    # print(idir)
    #
    # os.makedirs(fdir, exist_ok=True)
    #
    # if not os.path.isdir(idir):
    #     sys.exit("No input directory found.")
    #
    # for bc in barcodes.keys():
    #     print(f"Proceeding sample {barcodes[bc]}...")
    #     absite = {"A": 0,
    #               "C": 0,
    #               "G": 0,
    #               "T": 0,
    #               "N": 0}
    #
    #     idx = {"A": 0,
    #            "C": 1,
    #            "G": 2,
    #            "T": 3,
    #            "N": 4}
    #
    #     m = np.array(np.reshape(np.repeat(0, 5 ** 3), (5, 5, 5)))
    #
    #     fregex = f'({bc})' + '{e<=1}'
    #     rregex = str(Seq.reverse_complement(Seq(bc)))
    #     file = os.path.join(idir, (barcodes[bc]+"_final.fastq"))
    #
    #     for rec in SeqIO.parse(file, "fastq"):
    #
    #         indels = []
    #         f = re.search(fregex, str(rec.seq[0:9]))
    #         r = re.search(rregex, str(rec.seq[-10:-1]))
    #         if f and r:
    #             start = f.end()
    #             end = -(r.end()+2)
    #
    #             Fmm = rec.seq[55]
    #             As = rec.seq[66]
    #             Rmm = rec.seq[-55]
    #
    #
    #             rec = rec[start:end]
    #             #print(f"{rec.seq[55]} \t {rec.seq[-47]}")
    #             #print(f"{rec.seq[56:-48]}")
    #             if len(rec.seq[56:-47]) == 21:
    #                 #print(rec.seq[56:-47])
    #                 absite[str(rec.seq[66])] += 1
    #
    #                 Fmm = rec.seq[55]
    #                 As = rec.seq[66]
    #                 Rmm = rec.seq[-47]
    #
    #                 m[idx[Fmm], idx[As], idx[Rmm]] += 1
    #             else:
    #                 indels.append(rec)
    #     print(f"\tSaving indels seqs of sample {barcodes[bc]}...")
    #     SeqIO.write(indels, os.path.join(odir, f'{barcodes[bc]}_indels.fastq'), "fastq")
    #
    #     #print(f"In sample {barcodes[bc]} you lost {round(((ctr1-ctr2)/ctr1)*100, 4)}% of records")
    #     print(f"\tSaving heatplots of sample {barcodes[bc]}..")
    #     for nt in absite.keys():
    #         fig = plt.figure()
    #
    #         mm = (m[:, idx[nt], :] / (m[:, idx[nt], :]).sum())
    #
    #         df = pd.DataFrame(mm, columns=idx.keys(), index=idx.keys())
    #         fname = os.path.join(odir, f"{barcodes[bc]}_{nt}.csv")
    #         df.to_csv(path_or_buf=fname)
    #
    #         # hm = sns.heatmap(df, vmin = 0, vmax = 1, annot=True, cmap="YlGnBu" ).set_title(f"{barcodes[bc]}, {nt} (ratio is {round((absite[nt]/ctr3)*100, 2)})")
    #
    #         hm = sns.heatmap(df, vmin = 0, vmax = 1, annot=True, cmap="YlGnBu" ).set_title(f"{barcodes[bc]}, ({nt})")
    #
    #
    #         fig = hm.get_figure()
    #         fig_path = os.path.join(fdir, f"{barcodes[bc]}_{nt}.jpg")
    #         fig.savefig(fig_path)

