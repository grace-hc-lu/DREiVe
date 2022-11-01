package ComposeClusters_vert;
use Exporter;
use Cwd;
@ISA = ('Exporter');
@EXPORT = ('compose_clusters');

sub compose_clusters{
  my $properties = shift;
  my %properties = %$properties;
  my $win_length = $properties{"CLUSTER_LENGTH"};
  my $min_length = $properties{'MIN_PAT_LENGTH'};
  my $all_loci_pat = shift;
  my %all_loci_pat = %$all_loci_pat;
  my $patterns = shift;
  my %patterns = %$patterns;

  my %common_clusters = ();
  my %order_scores = ();
  my %borders = ();
  my %common_motifs = ();
  my %order_loci = ();

my %dist = qw(
equCab2  0.25
loxAfr3  0.3
canFam2  0.3
bosTau4  0.35
oryCun2  0.35
cavPor3  0.4
mm9      0.5
rn4      0.5
monDom5  0.8
taeGut1  1.1
galGal3  1.15
anoCar1  1.3
xenTro2  1.9
fr2      2.2
gasAcu1  2.1
oryLat2  2.4
danRer6  2.3

su       1.0
dasNov2  0.3
);

  $properties{'TIMEOUT'} = "";
  my $dir = cwd();
  my $motif_file = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."fasta.PC$properties{'MIN_PAT_LENGTH'}.pat1";  
  my $ps =  "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."ps";
  my $seq_number = 0;
  if ($properties{'SELECT'} eq 'all'){
    my $grep =  `/bin/grep '>' $motif_file`;
    my @grep = split(/\n/,$grep);
    $seq_number = scalar @grep;
  }elsif ($properties{'SELECT'} eq 'select'){
      $seq_number = $properties{'MINSPECIES'};
  }
  my $max_segment = (($properties{'REF_SEQ_END'} - $properties{'REF_SEQ_START'})/$properties{'SPLIT'}) + 1;
  for (my $segment = 1; $segment <= $max_segment; $segment ++){
    $motif_file = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."fasta.PC$properties{'MIN_PAT_LENGTH'}.pat"."$segment";  
    $cluster_file = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."fasta.PC$properties{'MIN_PAT_LENGTH'}.cluster$win_length".".$segment";
    $SIG{ALRM} = sub{die "timeout"};
    eval{
      alarm(3600);
      $pid = fork();
      unless ($pid){
        system "$properties{'DIR_ABS'}"."cgi-bin/Speedyclust/speedyclust_full $motif_file -J $seq_number -w $win_length -o $cluster_file < /dev/null > /dev/null 2>&1";
#        system "$properties{'DIR_ABS'}"."cgi-bin/Speedyclust/speedyclust_full $motif_file -J $seq_number -w $win_length -o $cluster_file";
#        system "$properties{'DIR_ABS'}"."cgi-bin/Promoclust/promoclust $motif_file -J $seq_number -w $win_length  -o $cluster_file";
        exit;
      }else{
        waitpid($pid,0);
      }
      alarm(0);
    };
    if ($@ =~ /timeout/){
      $properties{'TIMEOUT'} = "timeout";
      system "ps -eaf > $ps";
      $grep = `grep $pid $ps`;
      @grep = split (/\n/, $grep);
      foreach (@grep){
        /^\S+\s+(\d+)\s+(\d+)\s+/;
        $cpid = $1;
        $ppid = $2;
        kill ('INT',$cpid) if "$ppid" eq "$pid";
      }
      return(\%common_clusters,\%order_scores,\%borders,\%properties);
    }  
  }

  my @mdrosophilas = split(/,/,$properties{'MDROSOPHILAS'});
  my $dist_sum = 0; 
  my %dist_local = ();
  my $i = 1;
  foreach $specie (@mdrosophilas){
    $dist_sum += $dist{$specie};
    $dist_local{$i} = $dist{$specie};
    $i ++;
  }

  my $cutoff = 2*$properties{'MIN_PAT_LENGTH'};
  my $min_tokens = $properties{'MIN_PAT_LENGTH'};
  my $density_window = $properties{'DENSITY_WINDOW'};
  my $density_tokens = $properties{'DENSITY_TOKENS'};
  my $fasta = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."fasta";  
  my $fasta_upd = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."$win_length"."bp_$min_tokens"."bp_clusters.fasta.html";
  open (FASTAUPD,">$fasta_upd");
  my $tmp = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."tmp";  
  open(FASTA,"$fasta") || print "Can't open $fasta\n";
  $i = -1; my %fasta_seq = ();
  while(<FASTA>){
    chomp;
    if (/^>/){
      $i++;
    }else{
      $fasta_seq{$i} .= "$_";
    }
  }
  if ($properties{'SELECT'} eq 'all'){
    my $grep =  `/bin/grep -i ">" $fasta`;
    my @grep = split(/\n/,$grep);
    my $seq_number = scalar @grep;
  }elsif ($properties{'SELECT'} eq 'select'){
    my $seq_number = $properties{'MINSPECIES'};
  }

  my $f_motifs = 0; 
  my $f_coords = 0; 
  my $f_cutoff = 0;
  my $cluster_ID = 0;
  for ($segment = 1; $segment <= $max_segment; $segment ++){
    $cluster_file = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}$random"."fasta.PC$properties{'MIN_PAT_LENGTH'}.cluster$win_length".".$segment";  
    open (IN, "$cluster_file") || print "Can't open $cluster_file\n";
    CLUSTER:while (<IN>){
      chomp;
      if (/^MOTIF: (.*)$/){
        $motifs = $1;
        @motifs = split(/ <-> /,$motifs);
        $f_motifs = 1;
	%cluster_end = ();
      }elsif (/^{(.*)}$/){
        $coords = $1;
        @coords = split(/\)\(/,$coords); 
        next CLUSTER unless $coords[0] =~ /0\|(\d+),(\d+)/;
	foreach (@coords){
	  next unless /(\d+)\|(\d+),(\d+)/;
	  $specie = $1;
	  $start = $2;
	  $end = $3;
	  $cluster_end{$specie} = "$start-$end";
        }
        $f_coords = 1;
      }elsif (/^>(\d+)$/){
	$specie = $1;
      }elsif (/^(\(\d.*)$/ && $f_motifs && $f_coords){
	$motif_coord = $1;
        @motif_coord = split(/\)\(/,$motif_coord);
	$cluster_end{$specie} =~ /(\d+)-(\d+)/;
	$start = $1; $end = $2;
        foreach (@motif_coord){
	  next unless /(\d+),(\d+)/;
	  $motif = $1; 
          $coord = $2 + length($patterns{$motif});
	  $end = $coord if $end < $coord; 
	}
	$cluster_end{$specie} = "$start-$end";
	if ($specie == 0 && $end - $start < $cutoff){$f_motifs = 0; $f_coords = 0; next CLUSTER;}
	$f_cutoff = 1 if $specie == 0;
      }elsif (/\*\*\*/ && $f_motifs && $f_coords && $f_cutoff){
	$f_motifs = 0; $f_coords = 0; $f_cutoff = 0;
        $cluster_ID++;
	open (TMP,">$tmp"); %cluster_seq = ();
        for ($specie = 0; $specie <= @mdrosophilas; $specie++){
	  if ($cluster_end{$specie} =~ /(\d+)-(\d+)/){
	    $start = $1; $end = $2;
            if ($end <= length($fasta_seq{$specie})){
	      $cluster_seq = substr($fasta_seq{$specie},$start,$end-$start);
  	      $specie_name = $properties{'REF_SPECIES'} if $specie == 0;
	      $specie_name = $mdrosophilas[$specie-1] unless $specie == 0;
	    }else{
	      $cluster_seq = "AAAAA##################################################AAAAA";
	    }
	  }else{
	    $cluster_seq = "TTTTT##################################################TTTTT";
          }
	  print TMP ">$specie\n$cluster_seq\n";
	  $cluster_seq{$specie} = "$cluster_seq";
	}

        system "$properties{'DIR_ABS'}"."cgi-bin/Splash/splash -P standalone -a regular -q dna -i -j $seq_number -k $density_tokens -l $min_tokens -w $density_window -v -u -x 1000000 $tmp < /dev/null > /dev/null 2>&1";   
        my $add_pattern = 0; my %patterns_loc = ();
        my %motif_scores = (); my %motif_coord = ();
        my $motif_file = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."tmp.pat";
        open(MOTIFS, "$motif_file")||"Can't open $motif_file\n";
        MOTIFS:while(<MOTIFS>){
          chomp;
          unless(/;/){
            ($patID, $patStat, $pattern) = split(/\s/,$_);
            ($num_seq, $num_occur, $pat_length, $pval) = split(/\,/,$patStat);
            $num_seq =~ s/\[(\d+)/$1/;
            next if $num_seq < $seq_number;
            $patID =~ s/\[(\d+)\]/$1/;
            $add_pattern = 1;
          }elsif($add_pattern){
            ($offsets, $loci) = split(/;/,$_);
            @loci = split(/\]\[/,$loci);
            unless ($loci[0] =~ /\[0,/){
	      $add_pattern = 0; next MOTIFS;
            }
            %dist_l = %dist_local;
            foreach (@loci){
              /(\d+),(\d+)/;
              $sID = $1;
              $start = $2;
	      $motif_coord{$patID}{$sID}{$start} = 1;
              $motif_scores{$patID} += $dist_l{$sID} unless $sID == 0;
              $dist_l{$sID} = 0;
            }
            $motif_scores{$patID} = $motif_scores{$patID}/$dist_sum;
            $patterns_loc{$patID} = $pattern;
            $add_pattern = 0;
          }
        }

        my $add_cluster = 0;
        foreach $sID (sort {$a <=> $b} keys %cluster_end){
	  $cluster_end{$sID} =~ /(\d+)-(\d+)/;
	  $cluster_start = $1; $cluster_end = $2;       
          if ($sID == 0){  

#generate sequence of cluster length with motifs substituted with stretches of "X"            
            $seq = "N" x ($cluster_end - $cluster_start);
            %seq = ();
            foreach my $motif (keys %motif_coord){
              $x = "X" x length($patterns_loc{$motif});
              foreach my $loci (keys %{$motif_coord{$motif}{0}}){
		$start = $loci + $cluster_start;
                $common_clusters{$cluster_ID}{0}{$patterns_loc{$motif}}{$start} = 1;
		$motif_start = $loci;
                $motif_length = length($patterns_loc{$motif});
                $seq =~ s/(\w{$motif_start})\w{0,$motif_length}(.*)/$1$x$2/;
		while ($motif_length){
		  $seq{$motif_start} = $motif_scores{$motif} if ($seq{$motif_start} < $motif_scores{$motif});
		  $motif_length--;
                  $motif_start++;
                }
              }
            }

#calculate total number of motif lengths and islands
	    $motif_lengths = 0;
	    @seq = split(//,$seq);
            foreach my $coord (keys %seq){
	      $motif_lengths += $seq{$coord} if exists $seq{$coord};
            }
            $motif_lengths =~ s/(\d+)\.\d+/$1/;
	    if ($motif_lengths >= $cutoff){
  	      $cluster_end--;
              $borders{$cluster_ID}{0} = "$cluster_start:$cluster_end\n";    
    	      $order_loci{$cluster_start}{$cluster_end} = $cluster_ID;
    	      foreach(@motifs){$common_motifs{$_} = 1;}
              $add_cluster = 1;
	      $order_scores{$cluster_ID} = $motif_lengths;

	      @cluster_seq = split(//,$cluster_seq{0});
	      $cluster_color = "";
              foreach(@seq){
	        $n = shift(@cluster_seq);
	        $n = "<FONT COLOR = red>$n</FONT>" if /X/;
	        $cluster_color.= "$n";
	      }
              print FASTAUPD "Cluster $cluster_ID<BR>\n$properties{'REF_SPECIES'}<BR>\n$cluster_color<BR>\n";
      	    }else{
              delete($common_clusters{$cluster_ID});
              $add_cluster = 0;
	    }

          }elsif($add_cluster){   
            $seq = "N" x ($cluster_end - $cluster_start);
            foreach my $motif (keys %motif_coord){
              $x = "X" x length($patterns_loc{$motif});
              foreach my $loci (keys %{$motif_coord{$motif}{$sID}}){      
                $start = $loci + $cluster_start;       
                $common_clusters{$cluster_ID}{$sID}{$patterns_loc{$motif}}{$start} = 1;          
		$motif_start = $loci;
                $motif_length = length($patterns_loc{$motif});
                $seq =~ s/(\w{$motif_start})\w{0,$motif_length}(.*)/$1$x$2/;
              }
            }
	    $cluster_end--;
            $borders{$cluster_ID}{$sID} .= "$cluster_start:$cluster_end\n";

	    @seq = split(//,$seq);
	    @cluster_seq = split(//,$cluster_seq{$sID});
	    $cluster_color = "";
            foreach(@seq){
	      $n = shift(@cluster_seq);
	      $n = "<FONT COLOR = red>$n</FONT>" if /X/;
	      $cluster_color.= "$n";
	    }
            print FASTAUPD "$mdrosophilas[$sID-1]<BR>\n$cluster_color<BR>\n";
          }
        }
      }#if (/\*\*\*/)
    }#while (<IN>)
  }

  #add clusters that consist of single pattern
  foreach $patID (keys %all_loci_pat){
    foreach $loci (keys %{$all_loci_pat{$patID}{0}}){
      next unless (exists $common_motifs{$patID});

      $motif_score = 0; %dist_l = %dist_local;
      foreach my $sID (keys %{$all_loci_pat{$patID}}){
        $motif_score += $dist_l{$sID} unless $sID == 0;
        $dist_l{$sID} = 0;
      }
      $motif_score = ($motif_score/$dist_sum)*length($patterns{$patID});
      $motif_score =~ s/(\d+)\.\d+/$1/;

      if($motif_score >= $cutoff){
        $cluster_ID++; 
        $common_clusters{$cluster_ID}{0}{$patterns{$patID}}{$loci} = 1;
        $order_scores{$cluster_ID} = $motif_score;
        $loci_end = $loci + length($patterns{$patID}) - 1; 
        $borders{$cluster_ID}{0} = "$loci:$loci_end\n";
        $order_loci{$loci}{$loci_end} = $cluster_ID;
        $cluster_seq = substr($fasta_seq{0},$loci,length($patterns{$patID}));
        $cluster_color = "<FONT COLOR = red>"."$cluster_seq"."</FONT>";
        print FASTAUPD "Cluster $cluster_ID<BR>\n$properties{'REF_SPECIES'}<BR>\n$cluster_color<BR>\n";
        foreach $sID (sort {$a <=> $b} keys %{$all_loci_pat{$patID}}){
	  if ($sID){
            foreach $dcluster_start (sort {$a <=> $b} keys %{$all_loci_pat{$patID}{$sID}}){
              $common_clusters{$cluster_ID}{$sID}{$patterns{$patID}}{$dcluster_start} = 1;
              $dcluster_end = $dcluster_start + length($patterns{$patID}) - 1;  
              $borders{$cluster_ID}{$sID} .= "$dcluster_start:$dcluster_end\n";
              $cluster_seq = substr($fasta_seq{$sID},$dcluster_start,length($patterns{$patID}));
	      $cluster_color = "<FONT COLOR = red>"."$cluster_seq"."</FONT>";
              print FASTAUPD "$mdrosophilas[$sID-1]<BR>\n$cluster_color<BR>\n";
	    }
	  }
        }
      }
    }
  }


  @remove_ID = ();
  foreach $cluster_start (sort {$a <=> $b} keys %order_loci){
    foreach $cluster_end (sort {$b <=> $a} keys %{$order_loci{$cluster_start}}){
      $cluster_ID = $order_loci{$cluster_start}{$cluster_end};
      foreach $cur_cluster_start (sort {$a <=> $b} keys %order_loci){
	next if $cur_cluster_start < $cluster_start;
        last if $cur_cluster_start > $cluster_end;
        CUR: foreach $cur_cluster_end (sort {$b <=> $a} keys %{$order_loci{$cur_cluster_start}}){
	  if (($cur_cluster_start > $cluster_start && $cur_cluster_end < $cluster_end) || ($cur_cluster_start >= $cluster_start && $cur_cluster_end < $cluster_end) || ($cur_cluster_start > $cluster_start && $cur_cluster_end <= $cluster_end)){
	    $cur_cluster_ID = $order_loci{$cur_cluster_start}{$cur_cluster_end};
	    next CUR if $cur_cluster_ID == $cluster_ID;
	    foreach $sID (sort {$a <=> $b} keys %{$borders{$cluster_ID}}){
	      next unless $sID;
	      $s_start = 0; $s_end = 0;
	      $borders{$cluster_ID}{$sID} =~ /(\d+):(\d+)/;
	      $s_start = $1; $s_end = $2;
	      $cur_s_start = 0; $cur_s_end = 0;    
	      $borders{$cur_cluster_ID}{$sID} =~ /(\d+):(\d+)/;
	      $cur_s_start = $1; $cur_s_end = $2;  
              next unless $s_start && $s_end && $cur_s_start && $cur_s_end;
              next CUR unless (($s_start < $cur_s_start && $s_end > $cur_s_end) || ($s_start <= $cur_s_start && $s_end > $cur_s_end) || ($s_start < $cur_s_start && $s_end >= $cur_s_end));
	    }
	    push(@remove_ID,$cur_cluster_ID);
	  }
        }
      }
    }
  }
  foreach $ID (@remove_ID){
    delete($common_clusters{$ID});
    delete($order_scores{$ID});
    delete($borders{$ID});
  }

  return(\%common_clusters,\%order_scores,\%borders,\%properties);
}
