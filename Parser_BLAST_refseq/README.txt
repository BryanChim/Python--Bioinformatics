Parser_BLAST_refseq.py

****************************************
This program runs nucleotide blast for a user-specified queryset against a user-specified database file, and finds potential homologs by comparing (alignment/query length) vs. a user-defined threshold.

Command line arguments are accepted via argparse, in this order:
- name of file containing desired refseq ids to pull from query fasta file (ie. mouse_refseq_id.txt)
- name of query fasta file 
- name of database fasta file 
- desired name for blast output file (ie. mouse_human_blastoutfmt6)
- desired name for filtered output file (ie. filtered_out_06|07|08.txt)
- threshold for filtering blast output (minimum alignment length/query length ratio -- must be between 0.0 and 1.0)

Example Command Line to run program (for threshold of 0.8):

python refseq_blast_parser.py mouse_refseq_id.txt mouse.rna.fna human.rna.fna blast_out filtered_out_08.txt .8


****************************************
Here are the blast (-version 2.2.28) commands used, taken directly from the script:

os.system("makeblastdb -in " + args.dbfna + " -dbtype nucl -out blast_db")

os.system("blastn -db blast_db -query extracted_qrefseqs.fna -outfmt \"6 qseqid sseqid qlen slen length bitscore evalue\" -out " + args.blastout + " -evalue 1e-05")

Notes:

* As indicated above, the user specifies the database fasta file (human.rna.fna in this case) and the name of the blast output file (ie. "blast_out")

* Certain parts of the command that involve intermediate steps are hardcoded to simplify the process. They are:

- "blast_db" ~ name of the database file outputted from the makeblastdb command

- "extracted_qrefseqs.fna" ~ name of the fasta file containing extracted query sequences

- "-outfmt \"6 ... \"" ~ this blastn option specifies a custom tab-delimited format that only outputs the relevant information for each hit (query id, subject id, query length, subject length, alignment length, bitscore and evalue)

* The quotient, alignment length/query length is calculated for each hit -- if the result exceeds the user-defined threshold, that hit gets printed out to the filtered output file

 