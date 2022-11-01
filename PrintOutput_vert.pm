package PrintOutput_vert;
use Exporter;
@ISA = ('Exporter');    
@EXPORT = ('print_clusters');

sub print_clusters{
  my $order_length = shift;
  my $properties = shift;
  my $common_clusters = shift;
  my $borders = shift;
  my %order_length = %$order_length;
  my %properties = %$properties;
  my %common_clusters = %$common_clusters;
  my %borders = %$borders;
#  my ($clusterID,$pattern,$loci);
  my @mdrosophilas = split(/,/,$properties{'MDROSOPHILAS'});
  my @mdrosophilas = ("$properties{'REF_SPECIES'}",@mdrosophilas);
  my %mdrosophilas = ();
  my $i = 0;
  foreach (@mdrosophilas){
    $mdrosophilas{$i} = "$_";
    $i++;
  } 

  $cluster_number = 0;
  my $clusters_file =  "/$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.html";
  my $fasta_file = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.fasta.html";
  open (CLUSTERS, ">$clusters_file");
  print CLUSTERS "<HTML><HEAD><TITLE>CLUSTERS</TITLE></HEAD><BODY>\n";
  print CLUSTERS "Clusters of conserved motifs for $properties{$properties{'REF_SPECIES'}} gene $properties{'GENE'} limiting coordinates $properties{'REF_CHROM'}:$properties{'REF_SEQ_START'}-$properties{'REF_SEQ_END'}.<BR>Cluster length $properties{'CLUSTER_LENGTH'} bp and min motif length $properties{'MIN_PAT_LENGTH'} bp.<BR>Pattern density: $properties{DENSITY_TOKENS} matching nucleotides within window of $properties{DENSITY_WINDOW} bp.<BR>\n";
  foreach $clusterID (sort {$order_length{$b}<=>$order_length{$a}} keys %order_length){
    $cluster_number++;
    print CLUSTERS "<P>Cluster $cluster_number.<BR>\n";
    print CLUSTERS "Cluster score: $order_length{$clusterID}<BR>\n";
    print CLUSTERS "<TABLE FRAME=VOID RULES=ALL>\n";
    print CLUSTERS "<TR ALIGN=CENTER><TD>PATTERN</TD>";
    MDROSOPHILAS:foreach (@mdrosophilas){
      $fasta_seq =  "$properties{'URL'}"."cgi-bin/fasta.pl?CLUSTER=$clusterID&SPECIE=$_&FASTA=$fasta_file";
      print CLUSTERS "<TD><A HREF='$fasta_seq'>$_</A></TD>";
    }
    print CLUSTERS "</TR>\n";
    %print = ();
    %print_order = ();
    PATTERN:foreach $pattern (keys %{$common_clusters{$clusterID}{0}}){
      @dloci = ();
      foreach $loci (sort {$a <=> $b} keys %{$common_clusters{$clusterID}{0}{$pattern}}){
        push(@dloci,$loci);
      }
      @{$print{$pattern}} = @dloci;
      push (@{$print_order{$dloci[0]}},"$pattern");
    }

    foreach $loci0 (sort {$a <=> $b} keys %print_order){
      foreach $pattern(@{$print_order{$loci0}}){
          print CLUSTERS "<TR><TD><FONT SIZE=-1>$pattern</FONT></TD><TD>";
          foreach $loci (sort {$a <=> $b} @{$print{$pattern}}){
	    $loci += $properties{REF_SEQ_START};
            print CLUSTERS "$loci ";
          }
          print CLUSTERS "</TD>";
          foreach $i (sort {$a <=> $b} keys %mdrosophilas){
	    next unless $i;
            @loci = (); $f1 = 0;
            foreach $loci (sort {$a <=> $b} keys %{$common_clusters{$clusterID}{$i}{$pattern}}){   
	      push(@loci,$loci);
            }      
            print CLUSTERS "<TD>@loci</TD>";
          }
          print CLUSTERS "</TR>\n";  
      }
    }
    $borders{$clusterID}{0} =~ /(\d+):(\d+)/;
    $cluster_start = $1 + $properties{REF_SEQ_START};
    $cluster_end = $2 + $properties{REF_SEQ_START};
    print CLUSTERS "</TABLE>$properties{$properties{'REF_SPECIES'}} coordinates for cluster $cluster_start:$cluster_end<BR>\n";
  }
  close (CLUSTERS);

  unless ($cluster_number){
    $properties{'PRINT'} .= "Sorry, there are no clusters with max cluster length of $properties{'CLUSTER_LENGTH'} bp and min motif length $properties{'MIN_PAT_LENGTH'} bp.<BR>";
  }

  %properties;
}

