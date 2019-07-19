from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#TODO: sequence alignment in python

def match_seq_to_gene_list(infile, gene_list, output_folder):
    #returns the sequences that match the gene list
    #gene list should be generated from phosphorylation-cancer data
    #the distance matrix needs to be generated with software such as ugene pro
    seq_list = SeqIO.index(infile, 'fasta')
    output = []

    for index, row in gene_list.iterrows():
        tmp_name = 'Hs_' + row['Kinase'] + '.kin_dom'

        try:
            seq_rec = seq_list[tmp_name]
            tmp_rec = SeqRecord(seq_rec.seq,
                                id=row['Kinase'],
                                name=row['Kinase'],
                                description='')
            output.append(tmp_rec)
        except BaseException as e:
            pass

    SeqIO.write(output, output_folder + '/matched_gene_sequences.fasta', 'fasta')