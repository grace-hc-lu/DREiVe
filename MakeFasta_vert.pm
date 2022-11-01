package MakeFasta_vert;
use Exporter;
@ISA = ('Exporter');    
@EXPORT = ('make_fasta');

sub make_fasta{
  my $properties = shift;
  my %properties = %$properties;

  %properties = hello(\%properties);

  unless ($properties{'TEMPLATE'} eq "seq" || $properties{'TEMPLATE'} eq "file"){

    %properties = gene_info(\%properties) if $properties{'TEMPLATE'} eq "undef";
    return %properties if $properties{'PRINT'} =~ /Sorry/;

    %properties = seq_retrive(\%properties);
    return %properties if $properties{'PRINT'} =~ /Sorry/;

    %properties = mult_drosophila_seq(\%properties) if $properties{'MDROSOPHILAS'};

  }else{

    %properties = seq_template(\%properties);
    return %properties if $properties{'PRINT'} =~ /Sorry/;

  }

  %properties = lc_seq_masking(\%properties);

  %properties = print_fasta(\%properties);

  %properties;
}

sub hello{
  my $properties = shift;
  my %properties = %$properties;
  my $genome_annotation = "$properties{'ANNOTATION'}"."refFlat.txt";
  my $grep_genes = `/bin/grep -i "$properties{'GENE'}" $genome_annotation`;
  %properties;
}

#find info about Dmel gene
sub gene_info{
  my $properties = shift;
  my %properties = %$properties;
  my ($grep_genes, @grep_genes, $gene_start, $gene_end, $gene_name, $gene_id);  
  $genome_annotation = "$properties{'ANNOTATION'}"."refFlat.txt";
  $grep_genes = `/bin/grep -i "$properties{'GENE'}" $genome_annotation`;
  @grep_genes = split(/\n/,$grep_genes);
  foreach (@grep_genes){
    if (/^(\S+)\t(\S+)\t(\S+)\t(.)\t(\d+)\t(\d+)\t\d+\t\d+\t\d+\t\S+\t\S+/){
      $gene_name = $1;
      $gene_id = $2;
      $gene_chrom = $3;
      $gene_strand = $4;
      $gene_start = $5; 
      $gene_end = $6; 
      if ($gene_name =~ /^$properties{'GENE'}$/i || $gene_id =~ /^$properties{'GENE'}$/i){
        $properties{'REF_GENE_NAME'} = $gene_name; 
        $properties{'REF_GENE_ID'} = $gene_id; 
        $properties{'REF_CHROM'} = $gene_chrom;
        $properties{'REF_STRAND'} = $gene_strand; 
	if (exists($properties{'REF_GENE_START'}) && exists($properties{'REF_GENE_END'})){
	  $properties{'REF_GENE_START'} = $gene_start if $gene_start < $properties{'REF_GENE_START'};
	  $properties{'REF_GENE_END'} = $gene_end if $gene_end > $properties{'REF_GENE_END'};
        }else{
	  $properties{'REF_GENE_START'} = $gene_start;
	  $properties{'REF_GENE_END'} = $gene_end; 
	}
      }
    }
  }
  if ($properties{'REF_GENE_START'}){
    $properties{'REF_GENE_START'}++;
    $properties{'PRINT'} .= "Gene $properties{'GENE'} is located on $properties{$properties{'REF_SPECIES'}} chromosome $properties{'REF_CHROM'}:$properties{'REF_GENE_START'} - $properties{'REF_GENE_END'}.<BR>";
  }else{
    $properties{'PRINT'} .= "Sorry, there is no gene $properties{'GENE'} in $properties{$properties{'REF_SPECIES'}} genome annotation<BR>";
  }
  %properties;
}

sub seq_retrive{
  my $properties = shift;
  my %properties = %$properties;
  $gene_start = $properties{REF_GENE_START}-1;
  $gene_end = "$properties{REF_GENE_END}";
  $strand = "$properties{REF_STRAND}"; 
  $neighbors = "$properties{LIMIT_WITH_NEIGHBORS}";
  $min_flank = "$properties{MIN_FLANK}";
  $max_flank = "$properties{MAX_FLANK}";
  $set_flank = "$properties{SET_FLANK}";
  $gene_region = "$properties{GENE_REGION}";
  my ($seq_start,$seq_end);

  if ($properties{'TEMPLATE'} eq "undef"){
#find closest genes  
    my $genes_file = "$properties{'ANNOTATION'}"."$properties{'REF_CHROM'}.annotation";
    my ($current_start,$current_end);
    my $closest_start = 1000000000;
    my $closest_end = 0;
    if ($neighbors eq "yes"){
      open (GENES, "$genes_file");
      while (<GENES>){
	/^(\S+)\t(NM_.*?)\t\S+\t.\t(\d+)\t(\d+)\t/;
        $gene_name = $1; $gene_id = $2;
        $current_start = $3; $current_end = $4;
        $closest_end = $current_end if ($current_end < $gene_start && $closest_end < $current_end);
	$closest_end = $gene_start if ($current_start < $gene_start && $current_end > $gene_start && "$gene_name" ne "$properties{'REF_GENE_NAME'}" && "$gene_id" ne "$properties{'REF_GENE_ID'}");
        $closest_start = $gene_end if ($current_start < $gene_end && $current_end > $gene_end && "$gene_name" ne "$properties{'REF_GENE_NAME'}" && "$gene_id" ne "$properties{'REF_GENE_ID'}");
        $closest_start = $current_start if ($current_start > $gene_end && $closest_start > $current_start);
      }
      $properties{'PRINT'} .= "Limiting neighbors: closest upstream gene is located at $closest_end, closest downstream gene is located at $closest_start.<BR>";
    }

#define limiting coordinates
    if ($neighbors eq "yes"){
      if (($gene_start - $closest_end) >= $min_flank && ($gene_start - $closest_end) <= $max_flank){
        $seq_start = $closest_end + 1;
      }elsif (($gene_start - $closest_end) < $min_flank){
        $seq_start = $gene_start - $min_flank;
      }elsif (($gene_start - $closest_end) > $max_flank){
        $seq_start = $gene_start - $max_flank;
      }
      if (($closest_start - $gene_end) >= $min_flank && ($closest_start - $gene_end) <= $max_flank){
        $seq_end = $closest_start - 1;
      }elsif (($closest_start - $gene_end) < $min_flank){
        $seq_end = $gene_end + $min_flank;
      }elsif (($closest_start - $gene_end) > $max_flank){
        $seq_end = $gene_end + $max_flank;
      }
    }else{
      $seq_start = $gene_start - $set_flank;
      $seq_end = $gene_end + $set_flank;
    }

    if ($gene_region eq "up"){
      if ($strand eq "-"){
        $seq_start = $gene_end + 1;
      }else{
        $seq_end = $gene_start - 1;
      }    
    }elsif ($gene_region eq "down"){
      if ($strand eq "-"){
        $seq_end = $gene_start - 1;
      }else{
        $seq_start = $gene_end + 1;
      }
    }elsif ($gene_region eq "coding"){
      $seq_start = $gene_start;
      $seq_end = $gene_end;
    }

    if (($seq_end-$seq_start) > 500000 && $gene_region eq "full"){
      $gene_region = "coding";
      $properties{GENE_REGION} = "coding";
      $seq_start = $gene_start;
      $seq_end =$gene_end;
      $properties{'PRINT'} .= "Required sequences length for $properties{$properties{'REF_SPECIES'}} exceed 500 kb limit. DREiVe will limit analysis to the region between transcription start and transcription termination sites. <BR>";
    }
   
    if (($seq_end-$seq_start) > 500000 && ($gene_region eq "coding" || $gene_region eq "up" || $gene_region eq "down")){
      $properties{'PRINT'} .= "Sorry, required sequences length for $properties{$properties{'REF_SPECIES'}} exceed 500 kb limit. Please, specify another genomic region for DREiVe analysis.<BR>";
      return %properties;
    }

  }elsif ($properties{'TEMPLATE'} eq "def"){
    if ($properties{'RANGE'} =~ /(\w+):(\d+)-(\d+)/){
      $properties{'REF_CHROM'} = $1;
      $properties{'REF_STRAND'} = "+";    
      $seq_start = $2;
      $seq_end = $3;   
    }else{
      $properties{'PRINT'} .= "Sorry, region of $properties{$properties{'REF_SPECIES'}} genome that you wish to analyze was specified in a wrong format<BR>";
      return %properties;      
    }
  }

#retrieve sequence
  my ($fasta_file,$seq);
  my $chrom_seq = "";
  if ($properties{'EXCLUDE_UTR'} eq "on"){
    $fasta_file = "$properties{'GENOME'}"."$properties{'REF_CHROM'}"."_nt.rm.fasta";
  }else{
    $fasta_file = "$properties{'GENOME'}"."$properties{'REF_CHROM'}"."_nc.rm.fasta";
  }
  open (FASTA, "$fasta_file");
  while (<FASTA>){
    chomp;
    $chrom_seq .= "$_" unless /^>/;
  }

  $seq = substr($chrom_seq,($seq_start - 1),($seq_end - $seq_start + 1));
  my $seq_length = length ($seq);
  my $ref_seq = "$properties{REF_SPECIES}"."_SEQ";
  $properties{$ref_seq} = $seq;
  $properties{REF_SEQ_START} = $seq_start;
  $properties{REF_SEQ_END} = $seq_end;
  $properties{'PRINT'} .= "Limiting coordinates for $properties{$properties{'REF_SPECIES'}} sequence that was sent for DREiVe analysis $properties{'REF_CHROM'}: $seq_start-$seq_end, sequence length is $seq_length bp.<BR>";
  %properties;
}

sub mult_drosophila_seq{
  my $properties = shift;
  my %properties = %$properties;
  my @mdrosophilas = split(/,/,$properties{'MDROSOPHILAS'});
  my @min_tokens = split(/:/,$properties{'MIN_TOKENS'});
  my @min_tokens_sorted = sort {$a <=> $b} (@min_tokens);
  my $min_tokens = $min_tokens_sorted[0];
  foreach my $specie (@mdrosophilas){
    my $query_seq = "";
    my $seq_key = "$specie"."_SEQ";
    my $net1 = "$properties{'ALIGNMENT'}"."$specie/$properties{'REF_CHROM'}.$properties{REF_SPECIES}.$specie.fill"|| print "Can't open $net1\n";
    my $net2 = "$properties{'ALIGNMENT'}"."$specie/$properties{'REF_CHROM'}.$properties{REF_SPECIES}.$specie.gap"|| print "Can't open $net1\n";
    open (NET1, "$net1");
    while (<NET1>){
      /^(\d+) (\d+) (.*) (.) (\d+) (\d+)/;
      $ref_start = $0; $ref_end = 0; $query_chrom = ""; $query_strand = ""; $query_start = 0; $query_end = 0;
      $ref_start = $1;
      $ref_end = $ref_start + $2 - 1;
      $query_chrom = $3;
      $query_strand = $4;
      $query_start = $5;
      $query_end = $query_start + $6 - 1;
      next unless (($ref_start >= $properties{'REF_SEQ_START'} && $ref_start <= $properties{'REF_SEQ_END'})||($ref_end >= $properties{'REF_SEQ_START'} && $ref_end <= $properties{'REF_SEQ_END'})||($properties{'REF_SEQ_START'}>=$ref_start && $properties{'REF_SEQ_END'}<=$ref_end));
#print "$specie:$_     $ref_start $ref_end $query_chrom $query_strand $query_start $query_end\n";
      if ($ref_start < $properties{'REF_SEQ_START'}){
        open (NET2, "$net2");
        while (<NET2>){
	  /^(\d+) (\d+) (.*) (.) (\d+) (\d+)/;
	  $axt_ref_start = $1;
  	  $axt_ref_end = $axt_ref_start + $2;
          $axt_query_chrom = $3;
          $axt_query_strand = $4;
          if ("$query_chrom" eq "$axt_query_chrom" && "$query_strand" eq "$axt_query_strand"){
            $axt_query_start = $5;
            $axt_query_end = $axt_query_start + $6;
            if ($axt_ref_start <= $properties{'REF_SEQ_START'} &&  $properties{'REF_SEQ_START'} <= $axt_ref_end){
  	      if ($axt_query_strand eq "+" && $axt_query_start >= $query_start && $axt_query_start <= $query_end){
                $query_start = $axt_query_start;
		last;
	      }elsif ($axt_query_strand eq "-" && $axt_query_end >= $query_start && $axt_query_end <= $query_end){
                $query_end = $axt_query_end;
		last;
              }
	    }elsif ($axt_ref_start > $properties{'REF_SEQ_START'} && $axt_ref_start <= $properties{'REF_SEQ_END'}){
              $extra = $axt_ref_start - $properties{'REF_SEQ_START'};
	      if ($axt_query_strand eq "+" && $axt_query_start >= $query_start && $axt_query_start <= $query_end){
	        $query_start = $axt_query_start - $extra;
		last;
	      }elsif ($axt_query_strand eq "-" && $axt_query_end >= $query_start && $axt_query_end <= $query_end){
	        $query_end = $axt_query_end + $extra;
	        last;
              }
	    }
          }
        }
      }
#print "$specie:1:$query_start - $query_end\n";

      if ($ref_end > $properties{'REF_SEQ_END'}){
        open (NET2, "$net2");
        while (<NET2>){
	  /^(\d+) (\d+) (.*) (.) (\d+) (\d+)/;
	  $axt_ref_start = $1;
  	  $axt_ref_end = $axt_ref_start + $2;
          $axt_query_chrom = $3;
          $axt_query_strand = $4;
          if ("$query_chrom" eq "$axt_query_chrom" && "$query_strand" eq "$axt_query_strand"){
            $axt_query_start = $5;
            $axt_query_end = $axt_query_start + $6;
            if ($axt_ref_start <= $properties{'REF_SEQ_END'} &&  $properties{'REF_SEQ_END'} <= $axt_ref_end){
  	      if ($axt_query_strand eq "+" && $axt_query_end >= $query_start && $axt_query_end <= $query_end){
                $query_end = $axt_query_end;
		last;
	      }elsif ($axt_query_strand eq "-" && $axt_query_start >= $query_start && $axt_query_start <= $query_end){
                $query_start = $axt_query_start;
		last;
              }
	    }elsif ($axt_ref_start > $properties{'REF_SEQ_END'}){
              $extra = $axt_ref_start - $properties{'REF_SEQ_END'};
	      if ($axt_query_strand eq "+" && $axt_query_start >= $query_start && $axt_query_start <= $query_end){
	        $query_end = $axt_query_start - $extra;
		last;
	      }elsif ($axt_query_strand eq "-" && $axt_query_end >= $query_start && $axt_query_end <= $query_end){
	        $query_start = $axt_query_end + $extra;
	        last;
              }
	    }
	  }
        }
      }

#print "$specie:2:$query_start - $query_end\n";

      my $seq = "";
      if ($properties{'EXCLUDE_UTR'} eq "on"){
        $fasta_file = "$properties{'OTHER_GENOME'}"."$specie/$query_chrom"."_nt.rm.fasta";
      }else{
        $fasta_file = "$properties{'OTHER_GENOME'}"."$specie/$query_chrom"."_nc.rm.fasta";
      }
      open (FASTA, "$fasta_file");
      while (<FASTA>){
        $seq = $_ unless /^>/;
      }
      if ($query_start){
        $seq = substr($seq, $query_start - 1, $query_end - $query_start + 1);
      }else{
        $seq = substr($seq, 0, $query_end - $query_start + 1);
      }
      $seq = compliment_seq($seq) if $query_strand eq "-";
      $N = "N" x $min_tokens if $seq;
      $query_seq .= "$seq"."$N";
    }
    if (length($query_seq) > 500000){
      @ss = split(/,/,$properties{'MDROSOPHILAS'});
      $properties{'MDROSOPHILAS'} = "";
      foreach $ss(@ss){$properties{'MDROSOPHILAS'} .= "$ss," unless "$ss" eq "$specie";}
      $properties{'MINSPECIES'} --;
      $properties{'PRINT'} .= "Required sequences length for $specie exceed 500 kb limit. This sequence will be removed from your DREiVe search.<BR>\n";
      $properties{'PRINT'} .= "DREiVe will return clusters that are conserved in $properties{'MINSPECIES'} species.<BR>\n" if "$properties{SELECT}" eq "select";
    }elsif (length($query_seq) == 0){
      @ss = split(/,/,$properties{'MDROSOPHILAS'});
      $properties{'MDROSOPHILAS'} = "";
      foreach $ss(@ss){$properties{'MDROSOPHILAS'} .= "$ss," unless "$ss" eq "$specie";}
      $properties{'MINSPECIES'} --;
      $properties{'PRINT'} .= "There are no orthologous sequences for $specie. This specie will be removed from your DREiVe search.<BR>\n";
      $properties{'PRINT'} .= "DREiVe will return clusters that are conserved in $properties{'MINSPECIES'} species.<BR>\n" if "$properties{SELECT}" eq "select";
    }else{
      $properties{$seq_key} = $query_seq;
    }
  }
  %properties;
}

sub compliment_seq{
  my $seq = shift;
  my $compliment_seq = "";
  while ($seq){
    $nucl = chop($seq);
    if ($nucl eq "A"){
      $nucl = "T";
    }elsif ($nucl eq "T"){
      $nucl = "A";
    }elsif ($nucl eq "C"){
      $nucl = "G";
    }elsif ($nucl eq "G"){
      $nucl = "C";
    }
    $compliment_seq .= "$nucl";
  }
  return $compliment_seq;
}

sub seq_template{
  my $properties = shift;
  my %properties = %$properties;
  my $name = ""; 
  $properties{MDROSOPHILAS} = "";
  $properties{'REF_SPECIES'} = "";
  
  $fasta = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."fasta";
  open (FASTA,"$fasta");
  while (<FASTA>){
    if (/^>/){
      if (/^>((\w+)_SEQ)/){
        $name = $1;
        if ($properties{'REF_SPECIES'}){
          $properties{'MDROSOPHILAS'} .= "," if $properties{MDROSOPHILAS};
          $properties{'MDROSOPHILAS'} .= "$2";
        }else{
	  $properties{'REF_SPECIES'} = $2;
        }
      }else{
        $properties{PRINT} = "Sorry, you submitted sequences in a wrong format. Sequence name should have the following format: specie_SEQ (e.g. mm9_SEQ or rn4_SEQ).<BR>";   
        return %properties;   
      }
    }elsif (/[ATGCNatgcn#]+/){
      s/[^ATGCNatgcn#]//g;
      $length = length($_);
      if ($length > 500000){
        $properties{PRINT} = "Sorry, one of the sequences in you input file exceed 500 kb limit.";   
        return %properties;
      }
      $properties{$name} .= $_;
    }
  }
  $properties{REF_SEQ_START} = 1;
  my $ref_seq = "$properties{'REF_SPECIES'}"."_SEQ";
  $properties{REF_SEQ_END} = length($properties{$ref_seq});
  %properties;
}

#masking dinucleotide and mononucleotide repeats, N->#
sub lc_seq_masking{
  my $properties = shift;
  my %properties = %$properties;
  my @min_tokens = split(/:/,$properties{'MIN_TOKENS'});
  my @min_tokens_sorted = sort {$a <=> $b} (@min_tokens);
  my $min_tokens = $min_tokens_sorted[0];
  my @mdrosophilas = split(/,/,$properties{'MDROSOPHILAS'});
  @mdrosophilas = ("$properties{REF_SPECIES}", @mdrosophilas);
  foreach (@mdrosophilas){
    my $seq_key = "$_"."_SEQ";  
    my $seq = $properties{$seq_key}; 

$length = length($seq);
print "$seq_key - $length\n";

#mask Ns
    $seq =~ s/N/#/g;

#mask low compexity sequences
    my %mtri = (); my $tri = "";
    my %mdi = (); my $di = "";
    foreach my $first ("A","T","C","G"){
      foreach my $second ("A","T","C","G"){
        $di = "$first"."$second";
        foreach my $mfirst ("A","T","C","G"){
          $mdi{$di} .= ":$mfirst"."$second:" unless $mfirst eq $first;
        }
        foreach my $msecond ("A","T","C","G"){
          $mdi{$di} .= ":$first"."$msecond:" unless $msecond eq $second;
        }
        foreach my $third ("A","T","C","G"){
          $tri = "$first"."$second"."$third";
          foreach my $mfirst ("A","T","C","G"){
            $mtri{$tri} .= ":$mfirst"."$second"."$third:" unless $mfirst eq $first;
          }
          foreach my $msecond ("A","T","C","G"){
            $mtri{$tri} .= ":$first"."$msecond"."$third:" unless $msecond eq $second;
          }
          foreach my $mthird ("A","T","C","G"){
            $mtri{$tri} .= ":$first"."$second"."$mthird:" unless $mthird eq $third;
          }
        }
      }
    }

    $min_repeat = int ($min_tokens/2) + 1;
    $s = "#"x(($min_repeat-2)*2);
    $mseq = $seq;
    $coord = -2;
    while (length($mseq) > $min_repeat*2){
      if ($mseq =~ /^#/){
        $mseq =~ s/^(#*)(.*)/$2/;
        $coord += length($1);
      }else{
        $mseq =~ s/^(..)(.*)/$2/;
        $coord += 2;
        my $di = $1; my $subseq = $2;
        next if $di =~ /#/;
        my $i = 1; my $mutation = 0;
        while ($i < $min_repeat){
          $subseq =~ s/^(..)(.*)/$2/;
          my $next_di = $1;
          if ($next_di eq $di){
            $i++;
          }elsif ($mdi{$di} =~ /:$next_di:/){
            $mutation++;
            $i++;
          }else{
            last;
          }
        } 
        if ($mutation < ($min_repeat*2)/10 && $i == $min_repeat){
          $mask = substr($seq,$coord,$min_repeat*2);
          $mask =~ s/^(..).*(..)$/$1$s$2/;
          substr($seq,$coord,$min_repeat*2) = $mask;
        }
      }
    }

    my $min_repeat = int ($min_tokens/3) + 1;
    my $s = "#"x(($min_repeat-2)*3);
    my $mseq = $seq;
    my $coord = -3;
    while (length($mseq) > $min_repeat*3){
      if ($mseq =~ /^#/){
        $mseq =~ s/^(#*)(.*)/$2/;
        $coord += length($1);
      }else{
        $mseq =~ s/^(...)(.*)/$2/;
        $coord += 3;
        my $tri = $1; my $subseq = $2;
        next if $tri =~ /#/;
        my $i = 1; my $mutation = 0;
        while ($i < $min_repeat){
          $subseq =~ s/^(...)(.*)/$2/;
          my $next_tri = $1;
          if ($next_tri eq $tri){
            $i++;
          }elsif ($mtri{$tri} =~ /:$next_tri:/){
            $mutation++;
            $i++;
          }else{
            last;
          }
        }
        if ($mutation < ($min_repeat*3)/10 && $i == $min_repeat){
          $mask = substr($seq,$coord,$min_repeat*3);
          $mask =~ s/^(...).*(...)$/$1$s$2/;
          substr($seq,$coord,$min_repeat*3) = $mask;
        }
      }
    }

    for (my $i = 4; $i <= $min_tokens; $i++){
      my $min_repeat = int ($min_tokens/$i) + 1;
      $min_repeat = 3 if $min_repeat < 3;
      my $mseq = $seq;
      my $coord = 0 - $i;
      while (length($mseq) > $min_repeat*$i){
        if ($mseq =~ /^#/){
          $mseq =~ s/^(#*)(.*)/$2/;
          $coord += length($1);
        }else{
          $mseq =~ s/^((.){$i})(.*)/$3/;
          $coord += $i;
          my $tri = $1; my $subseq = $3;
          next if $tri =~ /#/;
          my $j = 1; 
          while ($j){
            $subseq =~ s/^((.){$i})(.*)/$3/;
            my $next_tri = $1;
            if ($next_tri eq $tri){
              $j++;
            }else{
              last;
            }
          }
          if ($j >= $min_repeat){
            my $mask = substr($seq,$coord,$j*$i);
            my $s = "#"x(($j-2)*$i);
            $mask =~ s/^((.){$i}).*((.){$i})$/$1$s$3/;
            substr($seq,$coord,$j*$i) = $mask;
	  }
        }
      }
    }

    $properties{$seq_key} = $seq;
  }#foreach (@mdrosophilas)
  %properties; 
}

sub print_fasta{
  my $properties = shift;
  my %properties = %$properties;
  my $mdrosophilas = "$properties{'REF_SPECIES'},"."$properties{MDROSOPHILAS}";
  my @mdrosophilas = split(/,/,$mdrosophilas);
  my @min_tokens = split(/:/,$properties{'MIN_TOKENS'});
  my @min_tokens_sorted = sort {$a <=> $b} (@min_tokens);
  my $min_tokens = $min_tokens_sorted[0];
  my $fasta = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."fasta";
  open (FASTA,">$fasta");
  foreach (@mdrosophilas){
    my $seq_key = "$_"."_SEQ";
    print FASTA ">$seq_key\n$properties{$seq_key}\n";
  }
  %properties;    
}
