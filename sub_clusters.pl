#!/usr/bin/perl

#require "Parse_Form.lib";
#&Parse_Form;
require "./cgi-lib.pl";
  my (%formdata,  # The form data
      %cgi_cfn,   # The uploaded file(s) client-provided name(s)
      %cgi_ct,    # The uploaded file(s) content-type(s).  These are
                  #   set by the user's browser and may be unreliable
      %cgi_sfn,   # The uploaded file(s) name(s) on the server (this machine)
      $ret,       # Return value of the ReadParse call.
      @data,
     );
  $cgi_lib::writefiles = "/tmp";
  $ret = &ReadParse(\%formdata,\%cgi_cfn,\%cgi_ct,\%cgi_sfn);

print "Content-type:text/html\n\n";
print "<HTML><HEAD><TITLE>Clusters</TITLE></HEAD><BODY>";
$start = $formdata{START};
$end = $formdata{END};
$file = $formdata{FILE};
@species = split(/\//,$formdata{SPECIES});


$f1 = 0;
$file =~ /(.*\/\d+bp_(\d+)bp_clusters).html/;
$seq_file = "$1".".seq.html";
$min_pat = $2;
open (SEQ,"$seq_file");
while (<SEQ>){
  $seq = $_ if $f1;
  last if $f1;
  $f1 = 1 if />$start - $end/;
}
print "$formdata{REF} sequence for selected conserved region (only motifs together with $min_pat bp of flanking sequences are shown):<BR>\n";
print "$seq<BR>\n";
foreach $species (@species){
  print "<A HREF = ''>$species</A> sequence<BR>";
}


print "<HR><P>Subset of clusters for selected conserved region:<BR>\n";
$print = ""; $f1 = 0;
open(CLUSTERS,"$file");
while (<CLUSTERS>){
  $f1 = 1 if />Cluster \d+\./;
  $print .= "$_" if $f1;
  if (/coordinates for cluster (\d+):(\d+)/){
    $cluster_start = $1;
    $cluster_end = $2;
    print "$print" if (($start <= $cluster_start && $cluster_start < $end) || ($start < $cluster_end && $cluster_end <= $end) || ($cluster_start <= $start &&  $cluster_end >= $end));
    $print = ""; $f1 = 0;
  }
}

print "</BODY></HTML>";vini


sub fasta_masked{
  local($enh_start, $enh_end, $cluster_file, $ref_start) = @_;
  local($fasta_file, $min_pat, $seq, $enh_seq, $motif, @motif_starts, $motif_start, $motif_length, $x, @enh_seq, $enh_masked, @seq, $n);
  $cluster_file =~ /(.*)\/\d+bp_(\d+)bp_clusters.html/;
  $fasta_file = "$1"."/fasta";
  $min_pat = $2;

  open (FASTA,"$fasta_file");
  while(<FASTA>){
    unless (/>/){
      $seq = $_;
      last;
    }
  }
  
  $enh_seq = substr($seq,($enh_start-$ref_start - $min_pat),($enh_end - $enh_start + 2*$min_pat + 1));
  $seq = "N" x ($enh_end - $enh_start + 2*$min_pat + 1);
  open (CLUSTERS,"$cluster_file");
  while (<CLUSTERS>){
    next unless /<TR><TD><FONT SIZE=-1>(.*?)<\/FONT><\/TD><TD>(.*?)<\/TD>/;
    $motif = $1;
    @motif_starts = split(/ /,$2);
    foreach $motif_start (@motif_starts){
      $motif_start =~ s/(\d+)/$1/;
      if ($motif_start >= $enh_start && $motif_start < $enh_end){
        $motif_length = length($motif) + $min_pat*2;
        $motif_start -= $enh_start;
        $x = "X" x $motif_length;
        $seq =~ s/(\w{$motif_start})\w{0,$motif_length}(.*)/$1$x$2/;
      }
    }
  }

  @enh_seq = split(//,$enh_seq);
  $enh_masked = "";
  @seq = split(//,$seq);
  foreach(@seq){
    $n = shift(@enh_seq);
    $n = "#" unless /X/;
    $enh_masked .= "$n";
  }

  return($enh_masked);
}
