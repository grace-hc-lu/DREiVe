#!/usr/bin/perl

use File::Copy;

@test_set = qw(GDF1);
$abs_dir = "/d/ls3/d/apache/dreive.cryst.bbk.ac.uk/";
$abs_html = "http://dreive.cryst.bbk.ac.uk/";

foreach $gene_name (@test_set){
  while (1){
    $code = int(rand 100000) + 1;
    $result = "$abs_dir"."htdocs/tmp/result"."$code";
    last unless (-e $result);
  }
  umask (000);
  mkdir("$result",0777);
  $features = $known{$gene_name};
  $gene_name =~ s/-/.ds./g;

  @param =   ("TEMPLATE-undef",
              "GENE-$gene_name",
              "LIMIT_WITH_NEIGHBORS-yes",
              "MIN_FLANK-5000",
              "MAX_FLANK-20000",
              "SET_FLANK-50000",
              "GENE_REGION-full",
              "EXCLUDE_UTR-on",
              "REF_SPECIES-hg19",
              "MDROSOPHILAS-equCab2,loxAfr3,canFam2,bosTau4,oryCun2,cavPor3,mm9,rn4,monDom5,taeGut1,galGal3,anoCar1,xenTro2",
              "SELECT-select",
              "MINSPECIES-10",
              "REPEAT_MASKER-on",
              "DENSITY_TOKENS-6",
              "DENSITY_WINDOW-8",
              "MIN_TOKENS-10:12",
              "CLUSTER_LENGTH-300:600:1000",
              "SCORE-2",
              "TOP-1.ds.10",
              "FEATURES-$features",
              "CODE-$code",
              "SPLIT-10000",
              "URL-$abs_html",
              "DIR-$abs_dir",
              "DIR_ABS-$abs_dir");

  print "$gene_name - results can be found at $result\n";
  system "main.vert.pl @param";

}
