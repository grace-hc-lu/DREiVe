#!/usr/bin/perl 

use Finished;
use Chrom;
use Contig;

$ref = "hg19";
@finished_species = (mm9);
@chrom_species = qw(equCab2 canFam2 bosTau4 oryCun2 rn4 monDom5 taeGut1 galGal3);
@contig_species = qw(loxAfr3 cavPor3 anoCar1 xenTro2);
@all = qw(hg19 equCab2 canFam2 bosTau4 oryCun2 mm9 rn4 monDom5 taeGut1 galGal3 loxAfr3 cavPor3 anoCar1 xenTro2);

foreach $species (@all){
  print "$species\n";
  $genome = "genomes/"."$species/";
  $annotation = "annotations/"."$species/";

#read annotation for each chromosome into %data and create annotation file for each chromosome
  %data = ();
  open (REFSEQ, "$annotation"."refFlat.txt");
  while (<REFSEQ>){
    $data{$1} .= $_ if /^\S+\tNM_.*?\t(\S+)\t.\t/;
  }
  foreach $chrom(keys %data){
    open (CHROM, ">$annotation"."$chrom.annotation");
    print CHROM "$data{$chrom}";
  }

#read sequence into %seq for each chromosome 
  opendir(DIR,"$genome");
  while(defined($file = readdir(DIR))){
    $chrom = "";
    %seq = ();
    next unless $file =~ /(.*)\.fa\.masked/;
    open (FASTA, "$genome"."$file");
    while (<FASTA>){
      chomp;
      $chrom = $1 if />(.*)/;
      $seq{$chrom} .= $_ unless /^>/;
    }

#mask nc and nt regions according to annotation files and creat nc and nt fasta files
    foreach $chrom (keys %seq){
      $nc_seq = $seq{$chrom}; $nt_seq = $seq{$chrom};
      $seq_length = length($seq{$chrom});
      open (REFSEQ, "$annotation"."$chrom.annotation");
      while (<REFSEQ>){
        if (/^\S+\tNM_.*?\t$chrom\t.\t\d+\t\d+\t(\d+)\t(\d+)\t\d+\t(\S+)\t(\S+)/){
          next if $2 > $seq_length;
          $CDS_start = $1;
          $CDS_end = $2;
          $exon_starts = $3;
          $exon_ends = $4;
          @exon_starts = split(/,/,$exon_starts);
          @exon_ends = split(/,/,$exon_ends);
          $i = 0;
          foreach $exon_start (@exon_starts){
	    next unless $exon_start =~ /\d/;
            $exon_end = $exon_ends[$i];
            $i++;
            $length = $exon_end - $exon_start;
            $N = "N" x $length;
            substr ($nt_seq, $exon_start, $length) = $N;      
            next if ($exon_start < $CDS_start && $exon_end < $CDS_start) || ($exon_start > $CDS_end && $exon_end > $CDS_end);
            $exon_start = $CDS_start if $exon_start < $CDS_start;
            $exon_end = $CDS_end if $exon_end > $CDS_end;
            $length = $exon_end - $exon_start;
            $N = "N" x $length;
            substr ($nc_seq, $exon_start, $length) = $N;
          }
        }
      }
      open (FASTA_NT, ">$genome"."$chrom"."_nt.rm.fasta");
      print FASTA_NT ">$chrom"."_nt\n$nt_seq\n";
      open (FASTA_NC, ">$genome"."$chrom"."_nc.rm.fasta");
      print FASTA_NC ">$chrom"."_nc\n$nc_seq\n";
    }
  }
}


#foreach $species (@all){
  $axt_dir = "alignments/$ref/$species/axtNet/";
  opendir(DIR,"$axt_dir");
  while(defined($file = readdir(DIR))){
    next if $file =~/short/;
    $axt_file = "$axt_dir"."$file";
    $axt_file_short = "$axt_file".".short";
    open (AXT, "$axt_file");
    open (AXT_SHORT, ">$axt_file_short");
    while (<AXT>){
      print AXT_SHORT if /^\d/; 
    }
  }
}


finished_alignment($ref, @finished_species);
chrom_alignment($ref, @chrom_species);
contig_alignment($ref, @contig_species);



