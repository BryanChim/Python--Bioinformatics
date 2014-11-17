Parser_GenBankFile.py accepts an input file in GenBank (.gbk) format that contains detailed information on a set of contigs, the genes contained within those contigs, start/end/orientation of each of those genes, and the organism.

the GBK file is parsed via an ordered series of regular expressions to extract and output to a file (EXAMPLE: hs_ref_filtered.txt):

<Organism>	<Chromosome>	<Contig>	<Gene_Name>	<Strand_Orientation>	<Start>	<End> 