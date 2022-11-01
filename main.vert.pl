#!/usr/bin/perl

use MakeFasta_vert;
use RunSplash_vert; 
use ProcessSplash_vert;
use ComposeClusters_vert;
use PrintOutput_vert;
use SP_vert;

foreach (@ARGV){
  ($key,$value) = split(/-/,$_);
  $value =~ s/\.ds\./-/g;
  $value =~ s/\.s\./ /g;
  $properties{$key} = $value;
}

$properties{'ANNOTATION'} = "annotations/"."$properties{'REF_SPECIES'}/";
$properties{'GENOME'} = "genomes/"."$properties{'REF_SPECIES'}/";
$properties{'ALIGNMENT'} = "alignments/"."$properties{'REF_SPECIES'}/";
$properties{'OTHER_GENOME'} = "genomes/";
$properties{'hg19'} = "human";
$properties{'equCab2'} = "horse";  
$properties{'loxAfr3'} = "elephant";  
$properties{'canFam2'} = "dog";  
$properties{'bosTau4'} = "cow";  
$properties{'oryCun2'} = "rabbit";  
$properties{'cavPor3'} = "guineapig";  
$properties{'mm9'} = "mouse";      
$properties{'rn4'} = "rat";      
$properties{'monDom5'} = "opossum";  
$properties{'taeGut1'} = "platypus";  
$properties{'galGal3'} = "chicken";  
$properties{'anoCar1'} = "lizard";  
$properties{'xenTro2'} = "frog";  
$properties{'fr2'} = "fugu";     
$properties{'gasAcu1'} = "stickleback";  
$properties{'oryLat2'} = "medaka";  
$properties{'danRer6'} = "zebrafish";  

$properties{'TITLE'} = "result"."$properties{'CODE'}/";

umask (000);
$prop_file = "$properties{'DIR'}"."htdocs/tmp/result"."$properties{'CODE'}"."/properties";
open (PROP,">$prop_file");
chmod (0666,$prop_file);
$old_fh = select(PROP);
$| = 1;
foreach $key (keys %properties){
  print "$key - $properties{$key}\n";
}
select($old_fh);

$tmp = "$properties{'DIR'}"."htdocs/tmp/result"."$properties{'CODE'}"."/main.html";
open (TMP,">$tmp");
chmod (0666,$tmp);
$old_fh = select(TMP);
$| = 1;
print "<HTML><HEAD><TITLE>Results</TITLE>";
print "<META HTTP-EQUIV = 'Refresh' CONTENT = '10'>";
print  "</HEAD><BODY><CENTER><H3>$properties{USER_TITLE}</H3></CENTER>This page will show the progress of your job and links to the result files.<BR>Your job can take from a few minutes to a few hours. <BR>E-mail will be sent to you once your job is completed and files with results will be saved on our server for 7 days.<BR>\n";
select($old_fh);

%properties = make_fasta(\%properties);
$fasta = "$properties{'URL'}"."tmp/$properties{'TITLE'}"."fasta";
print TMP "<P>$properties{'PRINT'}";
die if $properties{'PRINT'} =~ /Sorry/;
print TMP "<A HREF = $fasta>Sequences</A> in fasta format are prepared for SPLASH run.<BR>";
$properties{'PRINT'} = "";

%properties = run_splash(\%properties);
$splash = "$properties{'URL'}"."tmp/$properties{'TITLE'}"."fasta.pat";
if ($properties{'TIMEOUT'} =~ /timeout/){
  print TMP "<P>SPLASH run with min motif length of $properties{'MIN_PAT_LENGTH'} bp was not completed within 1 hour. You can try to increase this parameter and run DREiVe again.<P>";
  close (TMP);
}else{
  print TMP "<P>SPLASH motifs are generated with motif density of $properties{DENSITY_TOKENS} matching nucleotides within window of $properties{DENSITY_WINDOW} bp.<BR>";
}

open (PROP,">$prop_file");
chmod (0666,$prop_file);
$old_fh = select(PROP);
$| = 1;
foreach $key (keys %properties){
  print "$key - $properties{$key}\n";
}
select($old_fh);

my @cluster_length = split(/:/,$properties{"CLUSTER_LENGTH"});
my @cluster_length_sorted = sort {$a <=> $b} (@cluster_length);
my @min_pat_length = split(/:/,$properties{"MIN_TOKENS"});
my @min_pat_length_sorted = sort {$a <=> $b} (@min_pat_length);
print TMP "<P>Predicted conserved regions:<BR>";
MIN_PATTERN:foreach my $min_pat_length (@min_pat_length_sorted){
  $properties{"MIN_PAT_LENGTH"} = $min_pat_length;
  CLUSTER:foreach my $cluster_length (@cluster_length_sorted){
    $properties{"CLUSTER_LENGTH"} = $cluster_length;

    my ($all_loci_pat,$patterns) = process_motifs(\%properties);
    my %all_loci_pat = %$all_loci_pat; 
    my %patterns = %$patterns;
    $properties{'PRINT'} = "";

    my ($common_clusters,$order_scores,$borders,$properties) = compose_clusters(\%properties,\%all_loci_pat,\%patterns);
    my %common_clusters = %$common_clusters;
    my %order_scores = %$order_scores;
    my %borders = %$borders;
    %properties = %$properties;
    if ($properties{'TIMEOUT'} =~ /timeout/){
      print TMP "<P>Clustering motifs with min motif length of $properties{'MIN_PAT_LENGTH'} bp and maximal cluster length of $properties{'CLUSTER_LENGTH'} bp was not completed within 1 hour.<BR>";
      if ($min_pat_length == $min_pat_length[-1]){
        print TMP "You can try longer min motif length.<BR>" ;
      }else{
	print TMP "Longer min motif length will be tested.<BR>";
      }
      next CLUSTER;
    }

    %properties = print_clusters(\%order_scores,\%properties,\%common_clusters,\%borders);
    if ($properties{'PRINT'} =~ /Sorry/){
      print TMP "<P>$properties{PRINT}<BR>";
      $properties{'PRINT'} = "";  
      next;
    }else{   
      $clusters = "$properties{'URL'}"."tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.seq.html";
      print TMP "<P>cluster length $properties{'CLUSTER_LENGTH'} bp, min motif length $properties{'MIN_PAT_LENGTH'} bp<BR>";
      print TMP "<A HREF = $clusters>coordinates, sequences and putative binding sites</A> |";    
      %properties = scattered_plot(\%properties,\%borders,\%order_scores);
      $sp = "$properties{'URL'}"."tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.graph.html";
      print TMP " <A HREF = $sp>graph</A><BR>";   
    }
  }
}
close (TMP);

open (TMP,"$tmp");
$tmp_text = "";
while (<TMP>){
  s/Refresh//;
  $tmp_text .= "$_";
}
close (TMP);

open (TMP,">$tmp");
print TMP "$tmp_text";
print TMP "<P>Your job is completed!<BR>";
#if ($properties{GENE} =~ /\S/){
#  print TMP "<FORM ACTION = $properties{'URL'}"."cgi-bin/store_vert.pl METHOD = POST><BR>\n";
#  print TMP "<INPUT TYPE = submit VALUE = 'Store these results in DREiVe database'><BR>\n";
#  foreach $key (keys %properties){
#    print TMP "<INPUT TYPE = hidden NAME = $key VALUE = \"",$properties{$key},"\"><BR>\n" unless $key =~ /SEQ/;
#  }
#  print TMP "</FORM>";
#}
#print TMP "</BODY></HTML><BR>\n";
close (TMP);

$dir = "$properties{'DIR'}"."htdocs/tmp/result"."$properties{'CODE'}";
opendir(DIR,"$dir");
while(defined($file = readdir(DIR))){
  next unless ($file eq "fasta.pat" || $file =~ /fasta\.PC/);
  $file = "$dir/$file";
  unlink "$file";
}

if ($properties{'ADDRESS'}){
$name = 'DREiVe';
$my_address = 'a.sosinsky@mail.cryst.bbk.ac.uk';
$address = $properties{'ADDRESS'};
if ($properties{'USER_TITLE'}){
  $subject = "DREiVe results: $properties{'USER_TITLE'}";
}elsif($properties{'TEMPLATE'} eq 'undef'){
  $subject = "DREiVe results for gene $properties{'GENE'}";
}elsif($properties{'TEMPLATE'} eq 'def'){
  $subject = "DREiVe results for genomic region $properties{'RANGE'}";
}elsif($properties{'TEMPLATE'} eq 'file'){
  $subject = "DREiVe results";
}
$content = "Your job was completed and you can find results at $properties{'URL'}"."tmp/$properties{'TITLE'}"."main.html.";
open(SENDMAIL,"| /usr/lib/sendmail -t") || print TMP "Can not open SENDMAIL<BR>";
print SENDMAIL <<End_of_mail;
From: $name <$my_address>
To: $address
Subject: $subject
$content
End_of_mail
close(SENDMAIL);
}
