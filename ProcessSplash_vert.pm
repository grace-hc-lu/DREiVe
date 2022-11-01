package ProcessSplash_vert;
use Exporter;
@ISA = ('Exporter');
@EXPORT = ('process_motifs');

sub process_motifs{
  my $properties = shift;
  my %properties = %$properties;
  my($patID,$patStat,$pattern,$num_seq,$num_occur,$pat_length,$pval,$sID,$start,$offsets,$loci,@loci,$formated_file);
  my %local = ();
  my %all_loci_pat = ();
  my %patterns = ();
  my $add_pattern = 0;
  my @mdrosophilas = split(/,/,$properties{'MDROSOPHILAS'});
  my $seq_number = 0;
  if ($properties{'SELECT'} eq 'all'){
    $seq_number = scalar @mdrosophilas;
    $seq_number ++;
  }elsif ($properties{'SELECT'} eq 'select'){
    $seq_number = $properties{'MINSPECIES'};
  }

  my $motif_file = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."fasta.pat";
  open(MOTIFS, $motif_file);
  MOTIFS:while(<MOTIFS>){
    chomp;
    unless(/;/){
      ($patID, $patStat, $pattern) = split(/\s/,$_);
      ($num_seq, $num_occur, $pat_length, $pval) = split(/\,/,$patStat);
      $num_seq =~ s/\[(\d+)/$1/;
      next if $pat_length < $properties{'MIN_PAT_LENGTH'};
      next if $num_seq < $seq_number;
      $patID =~ s/\[(\d+)\]/$1/;
      $add_pattern = 1;
    }elsif($add_pattern){
      ($offsets, $loci) = split(/;/,$_);
      @loci = split(/\]\[/,$loci);
      unless ($loci[0] =~ /\[0,/){
	$add_pattern = 0; next MOTIFS;
      }
      foreach (@loci){
        /(\d+),(\d+)/;
        $sID = $1;
        $start = $2;
        $all_loci_pat{$patID}{$sID}{$start} = 1;
        $local{$sID}{$start}{$patID} = 1;
      }
      $patterns{$patID} = $pattern;
      $add_pattern = 0;
    }
  }

  my $segment;
  my %list_pat_ID = ();
  my $i = 0;
  for ($segment = 0; $segment <= $properties{'REF_SEQ_END'} - $properties{'REF_SEQ_START'}; $segment += $properties{'SPLIT'}){
    $i++;
    %list_pat_ID = ();
    $formated_file = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."fasta.PC$properties{'MIN_PAT_LENGTH'}.pat"."$i";
    open (OUT,">$formated_file");
    foreach $sID (sort {$a <=> $b} keys %local){
      %print_pat_ID = ();
      print OUT ">$sID\n";
      foreach $start (sort {$a <=> $b} keys %{$local{$sID}}){
        if ($sID eq "0" && $start >= $segment && $start <= $segment + $properties{'SPLIT'} + $properties{"CLUSTER_LENGTH"}){
	  $lpatID = ""; $lpattern = 0; @col_patterns = ();
          foreach $patID (sort {$a <=> $b} keys %{$local{$sID}{$start}}){
#print "$sID:$start - $patID - $patterns{$patID}\n";
	    if (length($patterns{$patID}) > $lpattern){
	      $lpattern = length($patterns{$patID});
	      $lpatID = $patID ;
	    }
	    push(@col_patterns, $patID);
          }

	  $lpatID = $list_pat_ID{$lpatID} if exists($list_pat_ID{$lpatID});
          print OUT "($lpatID,$start)";
#print "OUT1:($lpatID,$start)\n";
	  
	  foreach(@col_patterns){
	    unless (exists($list_pat_ID{$_})){
              $list_pat_ID{$_} = $lpatID;
#print "$_ - $list_pat_ID{$_}\n";
	    }else{
              print OUT "($list_pat_ID{$_},$start)" unless "$list_pat_ID{$_}" eq "$lpatID";
#print "OUT2:($list_pat_ID{$_},$start)\n" unless "$list_pat_ID{$_}" eq "$lpatID";
	    }
          }
        }elsif ($sID ne "0"){
          foreach $patID (sort {$a <=> $b} keys %{$local{$sID}{$start}}){
#print "$sID:$start - $patID\n";
	    if (exists $list_pat_ID{$patID}){
              print OUT "($list_pat_ID{$patID},$start)" unless exists $print_pat_ID{$list_pat_ID{$patID}}{$start};
#print "OUT$sID:$patID:($list_pat_ID{$patID},$start)\n"  unless exists $print_pat_ID{$list_pat_ID{$patID}}{$start};
              $print_pat_ID{$list_pat_ID{$patID}}{$start} = 1;
	    }
          }
        }
      }
      print OUT "\n";
    }
  }
  return(\%all_loci_pat,\%patterns);
}
