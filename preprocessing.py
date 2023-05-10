import pandas as pd
from Bio import SeqIO
import os
from configs import FileConfig, RunConfig
import regex
from functions import reverse_complement
from pprint import pprint


def preprocessing(runconfig: RunConfig, fileconfig: FileConfig):

    filename = fileconfig.file_path
    sample_name = filename.split('/')[1].split('_')[0]
    barcodes = runconfig.barcodes.dict()
    primers = runconfig.primers.dict()
    original_seq = runconfig.original_seq
    ref_seq = runconfig.ref_seq
    ref_seq_prefix = ref_seq[10:21]
    ref_seq_postfix = ref_seq[-16:-5]

    prefix_regex = f'({ref_seq_prefix})' + '{e<=1}'
    postfix_regex = f'({ref_seq_postfix})' + '{e<=1}'
    start = regex.search(prefix_regex, ref_seq).end()
    end = regex.search(postfix_regex, ref_seq).start()
    expected = len(ref_seq[start:end])

    odir = fileconfig.out_dir
    os.makedirs(odir, exist_ok=True)

    report = {}

    overall = 0
    all_strange = 0
    for bc in barcodes.keys():
        all, strange, matching, cropped, no_indel, indel = 0, 0, 0, 0, 0, 0
        outname = barcodes[bc]
        matching_length_seq = []
        indel_seq = []
        forward_primer_regex = f'({primers["forward"]})' + '{s<=1,e<=1,i<=1}'
        reverse_primer_regex = f'({primers["reverse"]})' + '{s<=1,e<=1,i<=1}'
        forward_barcode_regex = f'({bc})' + '{s<=1}'
        reverse_barcode_regex = f'({reverse_complement(bc)})' + '{s<=1}'
        original_seq_regex = f'({original_seq})' + '{s<=1,e<=1,i<=1}'
        original_seq_reverse_regex = f'({reverse_complement(original_seq)})' + '{s<=1,e<=1,i<=1}'
        prefix_regex = f'({ref_seq_prefix})' + '{s<=1,e<=1,i<=1}'
        postfix_regex = f'({ref_seq_postfix})' + '{s<=1,e<=1,i<=1}'
        print(f"Processing {barcodes[bc]} sample.")

        for rec in SeqIO.parse(filename, "fastq"):
            all += 1
            forward_matches = regex.findall(forward_barcode_regex, str(rec.seq[0:13]))
            reverse_matches = regex.findall(reverse_barcode_regex, str(rec.seq[-12::]))
            condition1 = len(forward_matches) > 0
            condition2 = len(reverse_matches) > 0
            if condition1 and condition2:
                matching += 1
                if bool(regex.search(forward_primer_regex, str(rec.seq))):
                    if not bool(regex.search(original_seq_regex, str(rec.seq))):
                        rec = rec
                    else:
                        continue
                elif bool(regex.search(reverse_primer_regex, str(rec.seq))):
                    if not bool(regex.search(original_seq_reverse_regex, str(rec.seq))):
                        rec = rec.reverse_complement(id="rc_"+rec.id, description="reverse complement")
                    else:
                        continue
                else:
                    strange += 1
                    continue

                if bool(regex.search(original_seq_regex, str(rec.seq))):
                    continue

                condition3 = bool(regex.search(prefix_regex, str(rec.seq)))
                condition4 = bool(regex.search(postfix_regex, str(rec.seq)))

                if condition3 and condition4:
                    cropped += 1
                    start = regex.search(prefix_regex, str(rec.seq)).end()
                    end = regex.search(postfix_regex, str(rec.seq)).start()

                    rec_final = rec[start:end]

                    if len(rec_final.seq) == expected:
                        no_indel += 1
                        matching_length_seq.append(rec_final)
                    else:
                        indel += 1
                        indel_seq.append(rec_final)

        SeqIO.write(matching_length_seq, os.path.join(odir, (outname + "_matching_length.fastq")), "fastq")
        SeqIO.write(matching_length_seq, os.path.join(odir, (outname + "_indels.fastq")), "fastq")

        overall += matching/all*100
        all_strange += strange/all*100

        report[barcodes[bc]] = {
            'number of matching barcodes': matching,
            '% of matching barcodes': round(matching/all*100, 2),
            'number of strange sequences': strange,
            '% of strange sequences': round(strange/all*100, 2),
            'number of cropped sequences': cropped,
            '% of cropped sequences in corresponding samples': round(cropped/matching*100, 2),
            'number of matching length reads': no_indel,
            '% matching length in cropped': round(no_indel/cropped*100, 2),
            'number of reads with indels': indel,
            '% indels in cropped': round(indel / cropped * 100, 2)
        }

    print(f'Preprocessed finished with {round(overall,2)}% matching and {round(all_strange, 2)}% '
          f'macthing but strange sequences.')
    print('Preprocessing finished.')

    pprint(report)

    report_path = os.path.join(odir, f'{sample_name}_report.csv')
    df = pd.DataFrame(report)
    df.to_csv(report_path)

        # SeqIO.write(selected_seqs, os.path.join(odir, (outname+".fastq")), "fastq")
        # print(f"Preprocessing step one is finished.")
        #
        # selected_seqs = []
        #
        # original = runconfig.original_seq
        # for rec in SeqIO.parse(os.path.join(odir, (outname+".fastq")), "fastq"):
        #     if original not in rec.seq:
        #         selected_seqs.append(rec)
        #
        # SeqIO.write(selected_seqs, os.path.join(odir, (outname + "_final.fastq")), "fastq")
        #
        # print(f"Preprocessing finished, there are {seq_ids} records in {barcodes[bc]} sample.")
