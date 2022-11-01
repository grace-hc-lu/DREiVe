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
print "<HTML><HEAD><TITLE>Potential binding transcription factors</TITLE></HEAD><BODY>";

$template = "$formdata{'FILE'}"."/patser.fasta";
$outfile = "$formdata{'FILE'}"."/patser.out";
open (SEQ, ">$template")||print "Can't write to file patser.fasta.";
print SEQ "seq \\$formdata{SEQ}\\\n";

%scores = (); %acc = ();
open (JASPAR, "./Jaspar/scores");
while (<JASPAR>){
  chmod();
  ($name,$score,$acc) = split(/\|/,$_);
  $scores{$name} = $score;
  $acc{$name} = $acc;
}
 
%TF = ();
foreach $name (keys %scores){
  next unless $name =~/\w/;
  local $ENV{"PATH"} = "./patser";
  $matrix = "./Jaspar/$name";
  $cut_off = 0.8*$scores{$name};
  $stderr = "/dev/null";
  `patser-v3d -m $matrix -w -f $template -c -ls $cut_off -A a:t 3 c:g 2 1>$outfile 2>$stderr`;

  open (PATSER, "$outfile");
  while (<PATSER>){
    if (/position=\s+(\d+)(\w?)\s+score=\s+\S+\s+ln\(p-value\)=\s+(\S+)\s+/){
      $loc = $1; $c = $2; $score = $3;
      if ($c){
	$strand = "+";
      }else{
	$strand = "-";
      }
      $loc += $formdata{START};
      push(@{$TF{$loc}},"$strand|$name|$score"); 
    }
  }
}

print "<H4><CENTER>List of putative binding sites for selected region.</CENTER></H4>"; 
print "Jaspar CORE vertebrate non-redundant collection of binding profiles with profile score threshold of 80% was used for scanning.<BR>";
print "<TABLE><TR ALIGN=CENTER><TD><B>Position</B></TD><TD><B>Strand</B></TD><TD><B>Factor name</B></TD><TD><B>ln(p-value)</B></TD>";
foreach $loc (sort {$a <=> $b} keys %TF){
  foreach (@{$TF{$loc}}){
    ($strand,$name,$score) = split(/\|/,$_);
    if ($acc{$name} =~ /,/ && $name =~ /::/){
      ($acc1,$acc2) = split(/,/,$acc{$name});
      ($name1,$name2) = split(/::/,$name);
      $html1 = "http://www.uniprot.org/uniprot/"."$acc1";
      $html2 = "http://www.uniprot.org/uniprot/"."$acc2";
      print "<TR ALIGN=CENTER><TD>$loc</TD><TD>$strand</TD><TD><A HREF = $html1>$name1</A>::<A HREF = $html2>$name2</A></TD><TD>$score</TD>";
    }elsif(defined($acc{$name})){
      $html = "http://www.uniprot.org/uniprot/"."$acc{$name}";
      print "<TR ALIGN=CENTER><TD>$loc</TD><TD>$strand</TD><TD><A HREF = $html>$name</A></TD><TD>$score</TD>";
    }else{
      print "<TR ALIGN=CENTER><TD>$loc</TD><TD>$strand</TD><TD>$name</TD><TD>$score</TD>";
    }
  }
  print "</TR>\n";
}


