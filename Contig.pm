package Contig;
use Exporter;
use Cwd;
@ISA = ('Exporter');
@EXPORT = ('contig_alignment');

sub contig_alignment{
my $ref = shift;
my @contig_species = @_;

foreach $species (@contig_species){
print "$species\n";
  $dir = "alignments/$ref/$species/";
  open (NET, "$dir/$ref\.$species.net")||print "Can not open net file\n";
  %gaps = (); %fills = (); $i = 0; $igap = 0; %fills_sorted = ();
  while(<NET>){
    $i++;
    if (/^net (.*) \d+/){
      $chrom = $1;
      %starts = (); 
    }elsif (/^(\s+)fill ((\d+) \d+ .* . \d+ \d+) /){
      $wsf = $1;
      $wsf = length($wsf);
      $start = $3;
      $fills{$wsf}{$chrom}{$start} = "$2\n";
      $starts{$wsf} = $start;
      if ($wsf == $wsg + 1 && $wsf > 1 && $i == $igap + 1){
        $fills{$wsf-2}{$chrom}{$starts{$wsf-2}} =~ /(\d+) (\d+) (.*) (.) (\d+) (\d+)/;
	delete($fills{$wsf-2}{$chrom}{$starts{$wsf-2}});
        $fill_start = $1;
	$fill_end = $1 + $2 - 1;
	$fill_info = $3;
	$fill_strand = $4;
        $fill_al_start = $5;   
        $fill_al_end = $5 + $6 - 1;   
        $last_gap = pop(@{$gaps{$chrom}});
        $last_gap =~ /(\d+) (\d+) .* . (\d+) (\d+)/;   
        $gap_start = $1;
        $gap_end = $1 + $2 - 1;
        $gap_al_start = $3;   
        $gap_al_end = $3 + $4 - 1;
        if ($fill_strand eq "+"){
          $length = $gap_start - $fill_start + 1;
          $length_al = $gap_al_start - $fill_al_start + 1;
          $fills{$wsf-2}{$chrom}{$fill_start} = "$fill_start $length $fill_info $fill_strand $fill_al_start $length_al\n";
          $length = $fill_end - $gap_end + 1;
          $length_al = $fill_al_end - $gap_al_end + 1;
          $fills{$wsf-2}{$chrom}{$gap_end} = "$gap_end $length $fill_info $fill_strand $gap_al_end $length_al\n";
        }elsif ($fill_strand eq "-"){
          $length = $gap_start - $fill_start + 1;
          $length_al = $fill_al_end - $gap_al_end + 1;
          $fills{$wsf-2}{$chrom}{$fill_start} = "$fill_start $length $fill_info $fill_strand $gap_al_end $length_al\n";
          $length = $fill_end - $gap_end + 1;
          $length_al = $gap_al_start - $fill_al_start + 1;
          $fills{$wsf-2}{$chrom}{$gap_end} = "$gap_end $length $fill_info $fill_strand $fill_al_start $length_al\n";
        }
        $starts{$wsf-2} = $gap_end;
      }
    }elsif (/^(\s+)gap (\d+ \d+ .* . \d+ \d+) /){
      $wsg = $1;
      $wsg = length($wsg);
      push (@{$gaps{$chrom}}, $2);
      $igap = $i;
    }
  }

  foreach $chrom (sort keys %{$fills{1}}){
    foreach $level (sort {$a <=> $b} keys %fills){
      foreach $start (sort {$a <=> $b} keys %{$fills{$level}{$chrom}}){ 
	  $fills_sorted{$chrom}{$start} = $fills{$level}{$chrom}{$start} if defined($fills{$level}{$chrom}{$start});
      }
    }

    open (FILLS, ">$dir/$chrom\.$ref\.$species.fill");
    open (GAPS, ">$dir/$chrom\.$ref\.$species.gap");
    foreach $gap(@{$gaps{$chrom}}){
      print GAPS "$gap\n";
    }  
    foreach $start1 (sort {$a <=> $b} keys %{$fills_sorted{$chrom}}){
      $fills_sorted{$chrom}{$start1} =~ /(\d+) (\d+) (.*) (.) (\d+) (\d+)/;
      $start2 = $1;
      $length = $2;
      $al_scaf = $3;
      $al_strand = $4;
      $al_start = $5;
      $al_length = $6; 
      if (($al_length - 100000) > $length){
        if ("$al_strand" eq "+"){
          $axt_file_short = "$dir/axtNet/$chrom\.$ref\.$species\.net.axt.short";
          open (AXT_SHORT,"$axt_file_short")||print "Can not open $axt_file_short file\n";; 
          %axt = (); @or_starts = (); @or_ends = (); @starts = (); @ends = ();
          while (<AXT_SHORT>){
	    /^\d+ .* (\d+) (\d+) (.*) (\d+) (\d+) (.)/;
	    $f1 = 0; $f2 = 0; 
            $axt_start = $1;
            $axt_end = $2;
	    $axt_al_scaf = $3;
            $axt_al_start = $4;
            $axt_al_end = $5;
            $axt_al_strand = $6;
            next unless ("$axt_al_scaf" eq "$al_scaf") && ("$al_strand" eq "$axt_al_strand");
	    $f1 = 1 if $axt_start >= $start2 &&  $axt_start <= $start2 + $length - 1;
	    $f1 = 1 if $axt_end >= $start2 &&  $axt_end <= $start2 + $length - 1;
    	    $f2 = 1 if $axt_al_start >= $al_start &&  $axt_al_start <= $al_start + $al_length - 1;
	    $f2 = 1 if $axt_al_end >= $al_start &&  $axt_al_end <= $al_start + $al_length - 1;
	    next unless $f1 && $f2;
	    $axt{$axt_al_start}  = "$_";
	    push(@starts, $axt_al_start);
	    push(@ends, $axt_al_end);
	    push(@or_starts, $axt_start);
	    push(@or_ends, $axt_end);
          }
	  $j = 0; $prev_gap = 0; $gap1 = 0; $gap2 = 0; 
          @or_starts = ($start2, @or_starts); 
          @or_ends = ($start2, @or_ends); 
          @starts = ($al_start, @starts); 
          @ends = ($al_start, @ends);
          $end2 = $start2 + $length - 1; 
          push(@or_starts, $end2); 
          push(@or_ends, $end2);
          $al_end = $al_start + $al_length - 1; 
          push(@starts, $al_end); 
          push(@ends, $al_end);
          foreach $start3 (@starts){
	    $gap = ($starts[$j+1] - $ends[$j]) - ($or_starts[$j+1] - $or_ends[$j]);;
            if ($gap > 100000 && $prev_gap < $gap){
	      $prev_gap = $gap;
      	      $gap1 = $starts[$j];
	      $gap2 = $starts[$j+1];
	    }
	    $j++;
	    last if $j == @starts - 1;
	  }
          if ($gap1 && $gap2){
	    $axt{$gap1} =~ /^\d+ .* \d+ (\d+) .* \d+ (\d+) /;
	    $length2 = $1 - $start2 + 1;
	    $al_length2 = $2 - $al_start + 1;
	    print FILLS "$start2 $length2 $al_scaf $al_strand $al_start $al_length2\n";

	    $axt{$gap2} =~ /^\d+ .* (\d+) \d+ .* (\d+) \d+ /;
	    $length2 = $start2 + $length - $1;
	    $al_length2 = $al_start + $al_length - $2;
	    print FILLS "$1 $length2 $al_scaf $al_strand $2 $al_length2\n";

          }elsif($length < 5000 && $al_length > 500000){

          }else{
            print FILLS "$fills_sorted{$chrom}{$start1}";
          }

        }elsif("$al_strand" eq "-"){
	  $chrom_seq = "genomes/$species/$al_scaf"."_nc\.rm\.fasta";
	  open (SEQ,$chrom_seq)||print "Can not open $chrom_seq file\n";
	  $seq = "";
	  while(<SEQ>){
	    chomp;
            $seq .= $_ unless /^>/;
          }
	  $chrom_length = length($seq);
	  $al_start = $chrom_length - ($al_start + $al_length);
          $axt_file_short = "$dir/axtNet/$chrom\.$ref\.$species\.net.axt.short";
          open (AXT_SHORT,"$axt_file_short")||print "Can not open $axt_file_short file\n";; 
          %axt = (); @or_starts = (); @or_ends = (); @starts = (); @ends = ();
          while (<AXT_SHORT>){
	    /^\d+ .* (\d+) (\d+) (.*) (\d+) (\d+) (.)/;
	    $f1 = 0; $f2 = 0; 
            $axt_start = $1;
            $axt_end = $2;
	    $axt_al_scaf = $3;
            $axt_al_start = $4;
            $axt_al_end = $5;
            $axt_al_strand = $6;
            next unless ("$axt_al_scaf" eq "$al_scaf") && ("$al_strand" eq "$axt_al_strand");
	    $f1 = 1 if $axt_start >= $start2 &&  $axt_start <= $start2 + $length - 1;
	    $f1 = 1 if $axt_end >= $start2 &&  $axt_end <= $start2 + $length - 1;
    	    $f2 = 1 if $axt_al_start >= $al_start &&  $axt_al_start <= $al_start + $al_length - 1;
	    $f2 = 1 if $axt_al_end >= $al_start &&  $axt_al_end <= $al_start + $al_length - 1;
	    next unless $f1 && $f2;
	    $axt{$axt_al_start}  = "$_";
	    push(@starts, $axt_al_start);
	    push(@ends, $axt_al_end);
	    push(@or_starts, $axt_start);
	    push(@or_ends, $axt_end);
          }
	  $j = 0; $prev_gap = 0; $gap1 = 0; $gap2 = 0; 
          @or_starts = ($start2, @or_starts); 
          @or_ends = ($start2, @or_ends); 
          @starts = ($al_start, @starts); 
          @ends = ($al_start, @ends);
          $end2 = $start2 + $length - 1; 
          push(@or_starts, $end2); 
          push(@or_ends, $end2);
          $al_end = $al_start + $al_length - 1; 
          push(@starts, $al_end); 
          push(@ends, $al_end);
          foreach $start3 (@starts){
	    $gap = ($starts[$j+1] - $ends[$j]) - ($or_starts[$j+1] - $or_ends[$j]);;
            if ($gap > 100000 && $prev_gap < $gap){
	      $prev_gap = $gap;
      	      $gap1 = $starts[$j];
	      $gap2 = $starts[$j+1];
	    }
	    $j++;
	    last if $j == @starts - 1;
	  }
          if ($gap1 && $gap2){
	    $axt{$gap1} =~ /^\d+ .* \d+ (\d+) .* \d+ (\d+) /;
	    $length2 = $1 - $start2 + 1;
	    $al_start2 = $chrom_length - $2;
	    $al_length2 = $2 - $al_start + 1;
	    print FILLS "$start2 $length2 $al_scaf $al_strand $al_start2 $al_length2\n";

	    $axt{$gap2} =~ /^\d+ .* (\d+) \d+ .* (\d+) \d+ /;
	    $length2 = $start2 + $length - $1;
	    $al_start2 = $chrom_length - ($al_start + $al_length);
	    $al_length2 = $al_start + $al_length - $2;
	    print FILLS "$1 $length2 $al_scaf $al_strand $al_start2 $al_length2\n";

          }elsif($length < 5000 && $al_length > 500000){

          }else{
            print FILLS "$fills_sorted{$chrom}{$start1}";
          }
        }
      }else{
        print FILLS "$fills_sorted{$chrom}{$start1}";
      }#if ($al_length - 100000) > $length
    }#foreach $start
  }#foreach $chrom
}#foreach $species
}
