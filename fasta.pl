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
print "<HTML><HEAD><TITLE>Sequence</TITLE></HEAD><BODY>";

$clusterID = $formdata{CLUSTER};
$specie = $formdata{SPECIE};
$f1 = 0; $f2 = 0;
open (FASTA,"$formdata{FASTA}");
while (<FASTA>){
  if ($f2 && $f1){
    $seq = $_ ;
    last;
  }
  if (/Cluster $clusterID\</){
    $f1 = 1;
  }elsif(/Cluster/){
    $f1 = 0;
  }
  $f2 = 1 if ($f1 && /$specie/);
}


print "$seq";
print "</BODY></HTML>";
