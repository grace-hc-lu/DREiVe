Create the following directories:
genomes/hg19/ and genomes/assembly_name/ for each specie to include
 list of species that can be included:
 equCab2 horse
 loxAfr3 elephant
 canFam2 dog
 bosTau4 cow
 oryCun2 rabbit
 cavPor3 guineapig
 mm9 mouse
 rn4 rat
 monDom5 opossum
 taeGut1 platypus
 galGal3 chicken
 anoCar1 lizard
 xenTro2 frog
		
annotations/hg19/ and annotations/assembly_name/ for each specie to include
alignments/hg19/assembly_name for each specie to include 	
htdocs/tmp
cgi-bin

Download genomic sequences, annotation and alignments from http://hgdownload.cse.ucsc.edu/downloads.html. 

For human genomic sequence download http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz into genomes/hg19. Then do the same for other species.

For human genome annotation download http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz into annotations/hg19. Then do the same for other species.

For alignments download .net files e.g. http://hgdownload.cse.ucsc.edu/goldenPath/hg19/vsMm9/hg19.mm9.net.gz into alignments/hg19/mm9 
Then download .axt files into alignments/hg19/mm9/axtNet/directory. DO the same for other species.

To decompress and unpack files 
gzip -d file.tar.gz
tar xvf file.tar

DREiVe scripts as well as Splash and Speedyclust can be downloaded from http://dreive.cryst.bbk.ac.uk/download/

Run format_genome.pl from DREiVe directory. This script transforms files into DREiVe-compatible formats.

List of files/directories to copy into cgi-bin:
stand_alone.vert.pl
main.vert.pl
MakeFasta_vert.pm
RunSpalsh_vert.pm
ProcessSplash_vert.pm
ComposeClusters_vert.pm
PrintOutput_vert.pm
SP_vert.pm
Splash/*
Speedyclust/*
patser/*
Jaspar/*
cgi-lib.pl
fasta.pl
full_length.pl
jaspar.pl
redraw_sp.pl
sub_clusters.pl

To generate graphics you also need to install Perl GD module from http://search.cpan.org/dist/GD/

Edit script stand_alone.vert.pl:
1. Write full path to DREiVe directory in $abs_dir and full URL in $abs_html. 
2. @test_set should contain gene names or RefSeq accession numbers for genes that you want to analyse.
3. Other parameters
LIMIT_WITH_NEIGHBORS - 'yes' or 'no' - analyzed region is limited by neighbouring genes
MIN_FLANK and MAX_FLANK - max and min limits for upstream and downstream sequences if LIMIT_WITH_NEIGHBORS is 'yes'
SET_FLANK - length of upstream and downstream sequences if LIMIT_WITH_NEIGHBORS is 'no'
GENE_REGION - 'up', 'down' or 'full'
MDROSOPHILAS - list of species included into analysis
MINSPECIES - min number of species that support each motif
DENSITY_TOKENS, DENSITY_WINDOW and MIN_TOKENS  - Splash parameter
CLUSTER_LENGTH - max length of cluster of short conserved motifs
SCORE - cut-off score




