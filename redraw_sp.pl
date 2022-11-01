#!/usr/bin/perl

#require "Parse_Form.lib";
#&Parse_Form;
use SP_vert;
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

@properties = split(/\n/,$formdata{properties});
foreach (@properties){
  /(\S+)=(\S+)/;
  $properties{$1} = $2;

}
$properties{REDRAW} = "yes";

@borders = split(/\n/,$formdata{borders});
foreach (@borders){
  /(\S+):(\S+)=(\S+)/;
  $borders{$1}{$2} = $3;
}

@order_length = split(/\n/,$formdata{order_length});
foreach (@order_length){
  /(\S+)=(\S+)/;
  $order_length{$1} = $2;
}

$properties{SCORE} = $formdata{SCORE};
$properties{LOC_DIR} = $formdata{LOC_DIR};
$properties{LOC_HTML} = $formdata{LOC_HTML};

unless ($properties{LOC_HTML}){
  $url = "$properties{URL}"."tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.graph.html";
}else{
  $url = "$properties{LOC_HTML}"."tmp.graph.html";
}
%properties = scattered_plot(\%properties,\%borders,\%order_length);

print "<HTML><HEAD><TITLE>Results</TITLE><META HTTP-EQUIV = 'Refresh' CONTENT = '0;URL=$url'></HEAD><BODY>";







