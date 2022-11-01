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
$prev_end = $formdata{START};
$last_end = $formdata{END};
$file = $formdata{FILE};
$cut_off = $formdata{SCORE};
$f1 = 0; $seq = "";

open(CLUSTERS,"$file");
while (<CLUSTERS>){
  last if $f1 && /full length masked sequences/;
  if (/\>\>(\d+) - (\d+)\</ && $f1){
    $start = $1;
    $end = $2;
  }elsif (/(.*)\<BR\>/ && $f1){
    $patterns = $1;
    if ($seq){
      $seq .= "#" x ($start - $prev_end - 1);
    }else{
      $seq .= "#" x ($start - $prev_end);
    }
    $prev_end = $end;
    $seq .= "$patterns";
  }
  $f1 = 1 if /sequence consrvation score above \<B\>$cut_off\<\/B\>/;
}
$seq .= "#" x ($last_end - $end);

print "$seq<BR>";
print "</BODY></HTML>";
