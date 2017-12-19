from Bio import SeqIO
import sys
import os

if not len(sys.argv) == 4:
    print "python\t" + sys.argv[0] + "\tinput_fasta\tnumber_of_seq_in_one_subfile\toutput_dir\n"
    sys.exit(0)

number = 0
subfile_seq_number = int(sys.argv[2])
file_number = 0
output_dir = sys.argv[3]
if not os.path.exists(output_dir):
    os.system('mkdir -p %s ' % output_dir)

file_name = output_dir + "/split_fa_" + str(file_number) + ".fa"
split_fa = open(file_name, "w")

with open(sys.argv[1], "r") as input_fasta:
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        number += 1
        seq_id = seq_record.id
        seq = seq_record.seq
        if number <= subfile_seq_number:
            split_fa.write(">{seq_id}\n{seq}\n".format(**locals()))
        else:
            split_fa.close()
            file_number += 1
            file_name = output_dir + "/split_fa_" + str(file_number) + ".fa"
            split_fa = open(file_name, "w")
            split_fa.write(">{seq_id}\n{seq}\n".format(**locals()))
            number = 1

split_fa.close()
