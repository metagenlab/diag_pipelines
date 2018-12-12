
def multiple_alignments2concatenated_alignments(fasta_files, out_name):

    from Bio import AlignIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Align import MultipleSeqAlignment
    # identification of all distinct fasta headers id (all unique taxons ids) in all fasta
    # storing records in all_seq_data (dico)
    taxons = []
    all_seq_data = {}
    for one_fasta in fasta_files:
        all_seq_data[one_fasta] = {}
        with open(one_fasta) as f:
            alignment = AlignIO.read(f, "fasta")
        for record in alignment:
            if record.id not in taxons:
                taxons.append(record.id)
            all_seq_data[one_fasta][record.id] = record


    # building dictionnary of the form: dico[one_fasta][one_taxon] = sequence
    concat_data = {}

    start_stop_list = []
    start = 0
    stop = 0

    for one_fasta in fasta_files:
        start = stop + 1
        stop = start + len(all_seq_data[one_fasta][list(all_seq_data[one_fasta].keys())[0]]) - 1
        start_stop_list.append([start, stop])
        for taxon in taxons:

            # check if the considered taxon is present in the record
            if taxon not in all_seq_data[one_fasta]:
                # if taxon absent, create SeqRecord object "-"*len(alignments): gap of the size of the alignment
                seq = Seq("-"*len(all_seq_data[one_fasta][list(all_seq_data[one_fasta].keys())[0]]))
                all_seq_data[one_fasta][taxon] = SeqRecord(seq, id=taxon)
            if taxon not in concat_data:
                concat_data[taxon] = all_seq_data[one_fasta][taxon]
            else:
                concat_data[taxon] += all_seq_data[one_fasta][taxon]

    # concatenating the alignments, writing to fasta file
    MSA = MultipleSeqAlignment([concat_data[i] for i in concat_data])
    with open(out_name, "w") as handle:
        AlignIO.write(MSA, handle, "fasta")


multiple_alignments2concatenated_alignments(snakemake.input, snakemake.output[0])
