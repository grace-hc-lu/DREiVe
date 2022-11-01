package SP_vert;
use Exporter;
use GD;
@ISA = ('Exporter');    
@EXPORT = ('scattered_plot');

sub scattered_plot{
  my $properties = shift;
  my %properties = %$properties;
  my $borders = shift;
  my %borders = %$borders;
  my $order_length = shift;
  my %order_length = %$order_length;
  my $cluster_length = $properties{'CLUSTER_LENGTH'};
  my $min_pat_length = $properties{'MIN_PAT_LENGTH'};
  my $cut_off = $properties{SCORE};
  unless ($properties{LOC_DIR}){
    $clusters_seq =  "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.seq.html";
    $clusters_file =  "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.html";
    $patser_file = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}";
    $sp = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.png";
    $graph = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.graph.html";
    $graph_link = "$properties{'URL'}"."/tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.png";
  }else{
    $clusters_seq =  "$properties{'LOC_DIR'}"."htdocs/tmp/$properties{'TITLE'}$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.seq.html";
    $clusters_file =  "$properties{'LOC_DIR'}"."$properties{'CLUSTER_LENGTH'}bp_$properties{'MIN_PAT_LENGTH'}bp_clusters.html";
    $patser_file = "$properties{'LOC_DIR'}"."htdocs/tmp/$properties{'TITLE'}";
    $sp = "$properties{LOC_DIR}"."tmp.png";
    $graph = "$properties{LOC_DIR}"."tmp.graph.html";
    $graph_link = "$properties{LOC_HTML}"."tmp.png";
  }

  $seq_key = "$properties{'REF_SPECIES'}"."_SEQ";
  my $dmel_seq = $properties{$seq_key};
  my @dmel_seq = split(//,$dmel_seq);
  my @mdrosophilas = split(/,/,$properties{'MDROSOPHILAS'});
  $species = "";
  foreach $mdrosophilas (@mdrosophilas){
    $species .= "$properties{$mdrosophilas}/";
  }
  @mdrosophilas = ("$properties{'REF_SPECIES'}",@mdrosophilas);
  $marg = 100;
  $scale = 10;
  $scale = 20 if length($dmel_seq) > 300000;
  $gc_window = 100;
  $gc_ruler = 100;
  $conservation_ruler = 100;
  $cluster_window = $cluster_length;
  $cluster_ruler = 100;

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

$dist_sum = 0;
foreach $specie (@mdrosophilas){
  $dist_sum += $dist{$specie};
}

#set image size and colors
  $widthpi = int(($properties{REF_SEQ_END} - $properties{REF_SEQ_START})/$scale);
  $imwidth = $widthpi+$marg*2;
  $imwidth = 1200 if $imwidth < 1000;
  $highpi = $gc_ruler+$conservation_ruler+$cluster_ruler+20;
  $im = GD::Image -> new ($imwidth,$highpi + $marg*2);
  $white = $im -> colorAllocate(255,255,255);
  $black = $im -> colorAllocate(0,0,0);#frame
  $grey = $im -> colorAllocate(200,200,200);#$grid
  $red = $im -> colorAllocate(255,49,24);#conservation
  $rose = $im -> colorAllocate(128,0,0);#features
  $TSS_blue = $im -> colorAllocate(0,107,132);#TSS&gene
  $coding_blue = $im -> colorAllocate(123,161,223);#coding shading
  $UTR_blue = $im -> colorAllocate(145,216,223);#UTR shading
  $blue = $im -> colorAllocate(198,236,215);#gap shading  
  $grey_blue = $im -> colorAllocate(220,220,240);#$repeats and low complexity


#set vertical ruler and grid
  $im -> line($marg,$marg,$marg,$gc_ruler+$conservation_ruler+$cluster_ruler+$marg,$black);
  $im -> line($marg+$widthpi,$marg,$marg+$widthpi,$gc_ruler+$conservation_ruler+$cluster_ruler+$marg,$black);
  $vert_grid = ((100*$scale)-($properties{REF_SEQ_START}%(100*$scale)))/$scale;
  $vert_coord = (int($properties{REF_SEQ_START}/(100*$scale))+1)*(100*$scale);
  while ($vert_grid <= $widthpi){
    $im -> line($vert_grid+$marg,$marg,$vert_grid+$marg,$gc_ruler+$conservation_ruler+$cluster_ruler+$marg,$grey);
    $vert_coord_kb = $vert_coord/1000;
    $im -> string(gdGiantFont,$vert_grid+$marg-20,$marg-15,"$vert_coord_kb",$black);
    $vert_grid += 100;
    $vert_coord += 100*$scale;
  }
  if ($properties{TEMPLATE} eq 'undef' || $properties{TEMPLATE} eq 'def'){
    $im -> string(gdGiantFont,$marg,30,"$properties{$properties{'REF_SPECIES'}} chromosome "."$properties{REF_CHROM},kb",$black);
  }elsif ($properties{TEMPLATE} eq 'file'){
    $im -> string(gdGiantFont,$marg,30,"$properties{$properties{'REF_SPECIES'}} sequence "."$properties{REF_CHROM},kb",$black);
  }


# mark '#' for Dmel sequence
  $j = 1;
  foreach (@dmel_seq){
    $im -> line($marg+$j/$scale,$marg,$marg+$j/$scale,$gc_ruler+$conservation_ruler+$cluster_ruler+$marg-1,$grey_blue) if /#/;
    $j++;
  }


#GC content Dmel
  $im -> string(gdGiantFont,$marg-35,$marg-10,"100%",$black);
  $im -> string(gdGiantFont,10,$marg+$gc_ruler/2-10,"GC content",$black);
  $im -> string(gdGiantFont,$marg-20,$marg+$gc_ruler-10,"0%",$black);
  $j = $gc_window/2;
  $prev_gc_content = 0;
  while($j < @dmel_seq - $gc_window/2){
    $gc_content = 0;
    for($i = $j - $gc_window/2; $i < $j + $gc_window/2; $i++){
      $gc_content++ if $dmel_seq[$i] =~ /G|C/;
    }
    $gc_content /= $gc_window;
    $im -> line($marg+$j/$scale-1,$marg+$gc_ruler*(1-$prev_gc_content),$marg+$j/$scale,$marg+$gc_ruler*(1-$gc_content),$black);
    $j += $scale;
    $prev_gc_content = $gc_content;
  }

#genes Dmel
  %exons = (); 
  $annotation = "$properties{'ANNOTATION'}"."$properties{REF_CHROM}.annotation";
  open (ANN, "$annotation");
  while (<ANN>){
    if (/^(\S+)\tNM_.*?\t\S+\t.\t\d+\t\d+\t\d+\t\d+\t\d+\t(\S+)\t(\S+)/){
      $name = $1;
      $exon_starts = $2; $exon_ends = $3; 
      @exon_starts = split(/,/,$exon_starts);
      @exon_ends = split(/,/,$exon_ends);
      $i = 0;
      foreach $exon_start (@exon_starts){
	next unless $exon_start =~ /\d/;
        $exon_end = $exon_ends[$i];
        $i++;
        if (($properties{REF_SEQ_START}-$scale <= $exon_start && $exon_start <= $properties{REF_SEQ_END}+$scale)||($properties{REF_SEQ_START}-$scale <= $exon_end && $exon_end <= $properties{REF_SEQ_END}+$scale)||($exon_start <= $properties{REF_SEQ_START} && $exon_end >= $properties{REF_SEQ_END})){  
          push (@{$exons{$name}},"$exon_start:$exon_end");
          $exon_start = $properties{REF_SEQ_START} if $exon_start < $properties{REF_SEQ_START}; 
          $exon_end = $properties{REF_SEQ_END} if $exon_end > $properties{REF_SEQ_END};   
          $im -> filledRectangle($marg+($exon_start-$properties{REF_SEQ_START})/$scale,$marg,$marg+($exon_end-$properties{REF_SEQ_START})/$scale,$gc_ruler+$conservation_ruler+$cluster_ruler+$marg-1,$UTR_blue);    
        }
      }
    }
  }

  %TSS = (); 
  $name_position = $marg-50;
  foreach $name (keys %exons){
    $grep = `/bin/grep -i "$name" $annotation`;
    @grep = split(/\n/,$grep);
    @starts = (); @ends = ();
    @sort_starts = (); @sort_ends = ();
    foreach (@grep){
      if (/^$name\tNM_.*?\t\S+\t(.)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t\d+\t\S+\t\S+/){
        $strand = $1;
        push(@starts,$2); push(@ends,$3); 
        $CDS_start = $4; $CDS_end = $5;
        foreach (@{$exons{$name}}){
          ($exon_start,$exon_end) = split(/:/,$_);
          if (($CDS_start<=$exon_start && $exon_start<=$CDS_end)||($CDS_start<=$exon_end && $exon_end<=$CDS_end)||($exon_start<$CDS_start && $exon_end>$CDS_end)){
            $exon_start = $CDS_start if $exon_start < $CDS_start; 
            $exon_start = $properties{REF_SEQ_START} if $exon_start < $properties{REF_SEQ_START}; 
            $exon_end = $CDS_end if $exon_end > $CDS_end;
            $exon_end = $properties{REF_SEQ_END} if $exon_end > $properties{REF_SEQ_END};     
            $im -> filledRectangle($marg+($exon_start-$properties{REF_SEQ_START})/$scale,$marg,$marg+($exon_end-$properties{REF_SEQ_START})/$scale,$gc_ruler+$conservation_ruler+$cluster_ruler+$marg-1,$coding_blue);  
	  }
        }                  
      }
    }
    @sort_starts = sort {$a <=> $b} @starts;
    @sort_ends = sort {$a <=> $b} @ends;
    $gene_start = $sort_starts[0];
    $gene_end = $sort_ends[-1];
    $gene_start_pi = ($gene_start - $properties{REF_SEQ_START})/$scale;
    $gene_start_pi = 0 if $gene_start_pi < 0;
    $gene_end_pi = ($gene_end - $properties{REF_SEQ_START})/$scale;
    $gene_end_pi = $widthpi if $gene_end_pi > $widthpi;
    $im -> filledRectangle($marg+$gene_start_pi,$marg-1,$marg+$gene_end_pi,$marg,$TSS_blue);

#mark transcription start site
    if ($strand eq "+"){
      $DTSS = ($gene_start - $properties{REF_SEQ_START})/$scale;
      next if defined $TSS{$name}{$DTSS};
      $name_position_hor = $marg+$DTSS+15;
      $name_position_hor = $marg if $name_position_hor < 0;
      $name_position_hor = $widthpi + $marg if $name_position_hor > $widthpi;
      $im -> string(gdGiantFont,$name_position_hor,$name_position,$name,$black);
      next if $DTSS < -1 || $DTSS > $widthpi+1;
      $im -> line($marg+$DTSS,$marg-40,$marg+$DTSS,$marg,$TSS_blue);
      $im -> line($marg+$DTSS,$marg-40,$marg+$DTSS+10,$marg-30,$TSS_blue);
      $im -> line($marg+$DTSS,$marg-20,$marg+$DTSS+10,$marg-30,$TSS_blue);
    }elsif ($strand eq "-"){
      $DTSS = ($gene_end - $properties{REF_SEQ_START})/$scale;
      next if defined $TSS{$name}{$DTSS};
      $name_position_hor = $marg+$DTSS-50;
      $name_position_hor = $marg if $name_position_hor < 0;
      $name_position_hor = $widthpi + $marg if $name_position_hor > $widthpi;
      $im -> string(gdGiantFont,$name_position_hor,$name_position,$name,$black);
      next if $DTSS < -1 || $DTSS > $widthpi+1;
      $im -> line($marg+$DTSS,$marg-40,$marg+$DTSS,$marg,$TSS_blue);
      $im -> line($marg+$DTSS,$marg-40,$marg+$DTSS-10,$marg-30,$TSS_blue);
      $im -> line($marg+$DTSS,$marg-20,$marg+$DTSS-10,$marg-30,$TSS_blue);
    }
    $name_position += 10;
    $name_position = $marg-50 if $name_position > $marg-30;     
    $TSS{$name}{$DTSS} = 1;
  }

#set horizontal ruler
  $im -> line($marg,$marg,$widthpi+$marg,$marg,$black);
  $im -> line($marg,$marg+$gc_ruler,$widthpi+$marg,$marg+$gc_ruler,$black);
  $im -> line($marg,$marg+$gc_ruler/2,$widthpi+$marg,$marg+$gc_ruler/2,$grey);
  $im -> line($marg,$marg+$gc_ruler+$conservation_ruler,$widthpi+$marg,$marg+$gc_ruler+$conservation_ruler,$black);
  $im -> line($marg,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler,$widthpi+$marg,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler,$black);  

#  $i = 0;
#  foreach ($properties{REF_SPECIE}){
#    $i++;
#    $im -> line($marg,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler+20*$i,$marg+$widthpi,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler+20*$i,$black);
#    $name_length = length($_);
#    $name_length = 10 if $name_length > 10;
#    $im -> string(gdGiantFont,$marg-10*$name_length,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler+20*$i-10,"$_",$black);
#  }

#set legend
  $im -> line(0,$highpi+$marg*1.5-10,0,$highpi+$marg*1.5+30,$TSS_blue);
  $im -> line(0,$highpi+$marg*1.5-10,10,$highpi+$marg*1.5,$TSS_blue);
  $im -> line(0,$highpi+$marg*1.5+10,10,$highpi+$marg*1.5,$TSS_blue);
  $im -> string(gdGiantFont,20,$highpi+$marg*1.5-8,"transcription",$black);
  $im -> string(gdGiantFont,20,$highpi+$marg*1.5+8,"start",$black);

  $im -> filledRectangle(150,$highpi+$marg*1.5+7,160,$highpi+$marg*1.5+13,$UTR_blue);
  $im -> string(gdGiantFont,170,$highpi+$marg*1.5-8,"untranslated",$black);
  $im -> string(gdGiantFont,170,$highpi+$marg*1.5+8,"regions",$black);

  $im -> filledRectangle(300,$highpi+$marg*1.5+7,310,$highpi+$marg*1.5+13,$coding_blue);
  $im -> string(gdGiantFont,320,$highpi+$marg*1.5-8,"coding",$black);
  $im -> string(gdGiantFont,320,$highpi+$marg*1.5+8,"sequences",$black);

  $im -> filledRectangle(450,$highpi+$marg*1.5+7,460,$highpi+$marg*1.5+13,$rose);
  $im -> string(gdGiantFont,470,$highpi+$marg*1.5,"features",$black);

  $im -> filledRectangle(600,$highpi+$marg*1.5+7,610,$highpi+$marg*1.5+13,$red);
  $im -> string(gdGiantFont,620,$highpi+$marg*1.5-8,"conserved",$black);
  $im -> string(gdGiantFont,620,$highpi+$marg*1.5+8,"regions",$black);

  $im -> filledRectangle(750,$highpi+$marg*1.5+7,760,$highpi+$marg*1.5+13,$grey_blue);
  $im -> string(gdGiantFont,770,$highpi+$marg*1.5-8,"interspersed repeats and",$black);
  $im -> string(gdGiantFont,770,$highpi+$marg*1.5+8,"low complexity sequences",$black);
  $im -> string(gdGiantFont,770,$highpi+$marg*1.5+24,"in D.melanogaster",$black);

#draw motifs
  $im -> string(gdGiantFont,10,$marg+$gc_ruler/2-10+$conservation_ruler,"motifs",$black);
  $im -> string(gdGiantFont,$marg-35,$marg+$gc_ruler+$conservation_ruler-13,"0"."bp",$black);

  $im -> string(gdGiantFont,5,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler/2-15,"conservation",$black);
  $im -> string(gdGiantFont,10,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler/2-5,"score",$black);
  $im -> string(gdGiantFont,$marg-20,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler-10,"0",$black);

  $max_length = 0;
  %motifs = ();
  %conservation_scores = ();
  open(CLUSTERS, "$clusters_file");
  while(<CLUSTERS>){
    chomp;    
    if (/^<TR><TD><FONT SIZE=-1>(.*?)<\/FONT><\/TD><TD>(.*?)(<\/TD>.*)/){
      $pattern = $1;
      $pattern_length_full = length($pattern);
      $starts = $2;
      $coords = $3;
      $pattern =~ s/\.//g;
      $pattern_length = length($pattern);
      $max_length = int($pattern_length*1.1) if $max_length < int($pattern_length*1.1);
      @coords = split(/<\/TD>/,$coords);
      $num = 0;
      $motif_score = 0;
      %dist_l = %dist;
      foreach (@coords){
        $specie = $mdrosophilas[$num];
	$motif_score += $dist_l{$specie} if /\d+/;
        $dist_l{$specie} = 0 if /\d+/;
	$num++;
      }
      $motif_score = $motif_score/$dist_sum;
      @starts = split(/ /,$starts);
      foreach $start (@starts){
	next unless $start =~ /\d+/;
        $start -= $properties{REF_SEQ_START};
        $motifs{$start} = $pattern_length if $pattern_length > $motifs{$start};
	$count = $pattern_length_full;
        while ($count){
          $conservation_scores{$start} = $motif_score if $conservation_scores{$start} < $motif_score;
          $start ++;
          $count --;
        }
      }
    }
  }
  $im -> string(gdGiantFont,$marg-35,$marg+$gc_ruler,"$max_length"."bp",$black);
  $im -> line($marg,$marg+$gc_ruler+$conservation_ruler*(1-$properties{'MIN_PAT_LENGTH'}/$max_length),$marg+$widthpi,$marg+$gc_ruler+$conservation_ruler*(1-$properties{'MIN_PAT_LENGTH'}/$max_length),$grey) if $max_length;
  foreach $start (sort {$a <=> $b} keys %motifs){
    $im -> line($marg+($start+$motifs{$start}/2)/$scale,$marg+$gc_ruler+$conservation_ruler*(1-$motifs{$start}/$max_length),$marg+($start+$motifs{$start}/2)/$scale,$marg+$gc_ruler+$conservation_ruler,$black);
  }

#############
#$png_im = $im -> png;
#$sp1 = "ku1."."$sp";
#open (IMAGE, "> $sp1")||print "Cannot open file $sp:$!";
#print IMAGE $png_im;
#close IMAGE;
###############

#calculate background correction
  $j = 0;
  $correction = 0;
  foreach (@dmel_seq){
    $correction += $conservation_scores{$j} if exists $conservation_scores{$j} && $_ ne "#";
    $j++;
  }
  $effective_length = $dmel_seq;
  $effective_length =~ s/#//g;
  $effective_length = length($effective_length);
  $correction /= $effective_length;

#draw cluster scores
  %conservation = ();  
  %plot_conservation = ();  
  $max_conservation = 0;
  $j = 0;
  while($j < @dmel_seq){
    $conservation = 0;
    if (substr($dmel_seq,$j,1) ne "#"){
      for($i = $j - $cluster_window/2; $i <= $j + $cluster_window/2; $i++){
        $conservation += $conservation_scores{$i}  if exists $conservation_scores{$i} && $i >= 0 && $i < @dmel_seq;
      }
      $conservation = ($conservation/$cluster_window)/$correction;
    }
    if ($j % $scale == 0){
      $plot_conservation{$j} = $conservation;
      $max_conservation = int($conservation*1.1) if $max_conservation < int($conservation*1.1);
    }
    push (@{$conservation{$conservation}},$j);
    $j++;
  }

  $im -> string(gdGiantFont,$marg-35,$marg+$gc_ruler+$cluster_ruler,"$max_conservation",$black);
  $prev_conservation = 0;
  foreach $j (sort {$a <=> $b} keys %plot_conservation){  
    $im -> line($marg+$j/$scale-1,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler*(1-$prev_conservation/$max_conservation),$marg+$j/$scale,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler*(1-$plot_conservation{$j}/$max_conservation),$black);
    $prev_conservation = $plot_conservation{$j};
  }
  $im -> line($marg,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler*(1-$cut_off/$max_conservation),$widthpi+$marg,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler*(1-$cut_off/$max_conservation),$grey);

####################
#$png_im = $im -> png;
#$sp2 = "ku2."."$sp";
#open (IMAGE, "> $sp2")||print "Cannot open file $sp:$!";
#print IMAGE $png_im;
#close IMAGE;
################

#select sequences with score above cut-off 
  $consrv_regions = "";
  @cut_off = (); 
  $inc = 0.5;
  @enh = (); 
  $cut_off_cur = $cut_off + 2*$inc;
  while($cut_off_cur >= $cut_off - 2*$inc){
    push (@cut_off, $cut_off_cur);
    $cut_off_cur -=  $inc;
  }

  foreach $cut_off_cur(@cut_off){
    @dmel_conservation = (); 
    for ($i = 0; $i < @dmel_seq; $i++){push (@dmel_conservation,"x");}  
    foreach $conservation (sort {$b <=> $a} keys %conservation){
      last unless $conservation;
      last if $conservation < $cut_off_cur;
      foreach $j (@{$conservation{$conservation}}){ 
        $dmel_conservation[$j] = "c";
      }
      $min_conservation = $conservation;
    } 

    $i = 0;
    $f1 = 1; $f2 = 0;
    @consrv_regions = ();
    $consrv_region_start = 0;
    $consrv_region_end = 0;
    foreach(@dmel_conservation){
      if (/c/ && $f1){
        $consrv_region_start = $i;
        $f1 = 0; $f2 = 1;
      }
      if ((/x/ && $f2) || (/c/ && $f2 && $i == @dmel_conservation - 1)){
        $consrv_region_end = $i;
        $consrv_region_end-- if /x/;
        $conserv_cluster_start = 0;
        $conserv_cluster_end = 0;
        $add_conserv_region = 0;
        foreach $cluster_ID (keys %borders){
    	  $borders{$cluster_ID}{0} =~ /(\d+):(\d+)/; 
          $border_start = $1;
          $border_end = $2;
          if (($border_start >= $consrv_region_start && $border_start <  $consrv_region_end) || ($border_end > $consrv_region_start && $border_end <= $consrv_region_end) || ($border_start <= $consrv_region_start && $border_end >= $consrv_region_end)){
	    unless ($conserv_cluster_start) {
              $conserv_cluster_start = $border_start; 
              $conserv_cluster_end = $border_end;
	    }else{
	      $conserv_cluster_start = $border_start if $conserv_cluster_start > $border_start;
              $conserv_cluster_end = $border_end if $conserv_cluster_end < $border_end;
	    }
	    $add_conserv_region ++;
          }
        }

        if ($add_conserv_region){
          $im -> filledRectangle($marg+$conserv_cluster_start/$scale,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler-1,$marg+$conserv_cluster_end/$scale,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler+1,$red) if $cut_off_cur == $cut_off;
          $conserv_cluster_start += $properties{'REF_SEQ_START'};   
          $conserv_cluster_end += $properties{'REF_SEQ_START'};      
          push(@consrv_regions,"$conserv_cluster_start-$conserv_cluster_end");
        }
        $f1 = 1; $f2 = 0;
      }
      $i++;
    }
    $prev_start = 0; $prev_end = 0; 
    if (@consrv_regions){
      $consrv_regions .= "\n<A NAME = '$cut_off_cur'><HR><P>$properties{$properties{REF_SPECIES}} sequences for clusters that overlap with the regions of sequence conservation score above <B>$cut_off_cur</B> (only conserved motifs together with $properties{'MIN_PAT_LENGTH'} bp flanking sequences are shown):<BR>\n";
      foreach (@consrv_regions){
        /(\d+)-(\d+)/;
        $start = $1;
        $end = $2;
        if ($start > $prev_end){
	  $print_start = $prev_start - $properties{'MIN_PAT_LENGTH'}; $print_end = $prev_end + $properties{'MIN_PAT_LENGTH'};
          if ($prev_end){
            $sub_clusters = "$properties{URL}"."cgi-bin/sub_clusters.pl?START=$print_start&END=$print_end&FILE=$clusters_file&REF=$properties{$properties{'REF_SPECIES'}}&SPECIES=$species";
            $consrv_regions .= ">$print_start - $print_end | <A HREF = $sub_clusters>subset of overlapping clusters</A> | ";
	    $seq = fasta_masked($prev_start,$prev_end,$clusters_file,$properties{'REF_SEQ_START'});
	    $seq_jaspar = $seq;
            $seq_jaspar =~ s/#/N/g;
            $jaspar = "$properties{URL}"."cgi-bin/jaspar.pl?START=$print_start&SEQ=$seq_jaspar&FILE=$patser_file";
            $consrv_regions .= "<A HREF = $jaspar>list of potential binding transcription factors</A><BR>\n";
	    $consrv_regions .= "$seq<BR>\n";
	  }
	  push(@enh,"$print_start-$print_end") if ($prev_end && $cut_off_cur == $cut_off);
          $prev_start = $start;
          $prev_end = $end;
        }else{
	  $prev_end = $end;
        }
      }
      $print_start = $prev_start - $properties{'MIN_PAT_LENGTH'}; $print_end = $prev_end + $properties{'MIN_PAT_LENGTH'};
      if ($prev_end){
        $sub_clusters = "$properties{URL}"."cgi-bin/sub_clusters.pl?START=$print_start&END=$print_end&FILE=$clusters_file&REF=$properties{$properties{'REF_SPECIES'}}&SPECIES=$species";
        $consrv_regions .= ">$print_start - $print_end | <A HREF = $sub_clusters>subset of overlapping clusters</A> | ";
        $seq = fasta_masked($prev_start,$prev_end,$clusters_file,$properties{'REF_SEQ_START'});
        $seq_jaspar = $seq;
        $seq_jaspar =~ s/#/N/g;
        $jaspar = "$properties{URL}"."cgi-bin/jaspar.pl?START=$print_start&SEQ=$seq_jaspar&FILE=$patser_file";
        $consrv_regions .= "<A HREF = $jaspar>list of potential binding transcription factors</A><BR>\n";
        $consrv_regions .= "$seq<BR>\n";
      }
      $full_length = "$properties{URL}"."cgi-bin/full_length.pl?FILE=$clusters_seq&SCORE=$cut_off_cur&START=$properties{'REF_SEQ_START'}&END=$properties{'REF_SEQ_END'}";
      $consrv_regions .= "<A HREF = $full_length>full length masked sequences</A><BR>\n";
      push(@enh,"$print_start-$print_end") if ($prev_end && $cut_off_cur == $cut_off);
    }
  }

####################
#$png_im = $im -> png;
#$sp3 = "ku3."."$sp";
#open (IMAGE, "> $sp3")||print "Cannot open file $sp:$!";
#print IMAGE $png_im;
#close IMAGE;
################

  unless ($properties{REDRAW} eq "yes"){
    open (CLUSTERS, ">$clusters_seq");
    if ($properties{'TEMPLATE'} eq "undef"){
      print CLUSTERS "Predicted regulatory regions for $properties{$properties{REF_SPECIES}} gene $properties{'GENE'} (limiting coordinates on $properties{'REF_CHROM'}:$properties{'REF_SEQ_START'}-$properties{'REF_SEQ_END'}).<BR>\n";
    }elsif ($properties{'TEMPLATE'} eq "def"){
      print CLUSTERS "Predicted regulatory regions for $properties{$properties{REF_SPECIES}} $properties{'REF_CHROM'}:$properties{'REF_SEQ_START'}-$properties{'REF_SEQ_END'}.<BR>\n";
    }elsif ($properties{'TEMPLATE'} eq "file"){
      print CLUSTERS "Predicted regulatory regions.<BR>\n";
    }

    print CLUSTERS "Cluster length $properties{'CLUSTER_LENGTH'} bp and min motif length $properties{'MIN_PAT_LENGTH'} bp.<BR>Motif density: $properties{DENSITY_TOKENS} matching nucleotides within window of $properties{DENSITY_WINDOW} bp.<BR>\n";
    print CLUSTERS "Sequence conservation score cut-off:";
    foreach $score(@cut_off){
      print CLUSTERS " <A HREF = '#$score'>$score</A> |";
    }
    print CLUSTERS "<BR>";
    print CLUSTERS "$consrv_regions";
  }



#draw top clusters
#  @top = split(/,/,$properties{TOP});
#  %top = ();
#  foreach (@top){
#    if (/(\d+)-(\d+)/){
#      $first = $1; $last = $2;
#      for ($i = $first; $i <= $last; $i++){$top{$i} = 1;}
#    }elsif(/(\d+)/){
#      $top{$1} = 1;
#    }
#  }
#  $count_clusters = 0;
#  ID:foreach $clusterID (sort {$order_length{$b}<=>$order_length{$a}} keys %order_length){
#    $count_clusters++;
#    next unless defined $top{$count_clusters};
#    @borders = split(/\n/,$borders{$clusterID}{'0'});
#    foreach (@borders){
#	/(\d+):(\d+)/;
#        $cluster_start = $1;
#        $cluster_end = $2;
#        $im -> line($marg+$cluster_start/$scale,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler+20*($specie+1),$marg+$cluster_end/$scale,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler+20*($specie+1),$white);
#        $im -> string(gdGiantFont,$marg+($cluster_start+($cluster_end-$cluster_start)/2-10)/$scale,$marg+$gc_ruler+$conservation_ruler+$cluster_ruler+20*($specie+1)-10,"$count_clusters",$black);
#    }
#  }


#draw features
  $features = $properties{'FEATURES'};
  @features = split(/:/,$features);
  foreach(@features){
    /(\d+)-(\d+)/;
    $feature_start = $1;
    $feature_end = $2;
    $im -> filledRectangle($marg+($feature_start-$properties{REF_SEQ_START})/$scale,$marg-1,$marg+($feature_end-$properties{REF_SEQ_START})/$scale,$marg+1,$rose);
  }


#print image
  $png_im = $im -> png;
  open (IMAGE, "> $sp")||print "Cannot open file $sp:$!";
  print IMAGE $png_im;
  close IMAGE;
  open (GRAPH, ">$graph")||print "Cannot open file $graph:$!";
  print GRAPH "<HTML><HEAD><TITLE>GRAPH</TITLE></HEAD><BODY>\n";
  if ($properties{'TEMPLATE'} eq "undef"){
    print GRAPH "Predicted regulatory regions for $properties{$properties{REF_SPECIES}} gene $properties{'GENE'} (limiting coordinates on $properties{'REF_CHROM'}:$properties{'REF_SEQ_START'}-$properties{'REF_SEQ_END'}).<BR>\n";
  }elsif ($properties{'TEMPLATE'} eq "def"){
    print GRAPH "Predicted regulatory regions for $properties{$properties{REF_SPECIES}} $properties{'REF_CHROM'}:$properties{'REF_SEQ_START'}-$properties{'REF_SEQ_END'}.<BR>\n";
  }elsif ($properties{'TEMPLATE'} eq "file"){
    print GRAPH "Predicted regulatory regions.<BR>\n";
  }

  print GRAPH "Cluster length $properties{'CLUSTER_LENGTH'} bp and min motif length $properties{'MIN_PAT_LENGTH'} bp.<BR>\n";
  print GRAPH "Motif density: $properties{DENSITY_TOKENS} matching nucleotides within window of $properties{DENSITY_WINDOW} bp.<BR>\n";

  print GRAPH "<FORM ACTION=$properties{'URL'}/cgi-bin/redraw_sp.pl METHOD=POST>\n";

  $properties = "";
  foreach $key (keys %properties){
    $properties .= "$key=$properties{$key}\n" if $key && $properties{$key};
  }
  print GRAPH "<INPUT TYPE = hidden NAME = properties VALUE = \"",$properties,"\">\n";

  $borders = "";
  foreach $cluster_ID (keys %borders){
    foreach $specie (keys %{$borders{$cluster_ID}}){
      $borders .= "$cluster_ID:$specie=$borders{$cluster_ID}{$specie}\n";
    }
  }
  print GRAPH "<INPUT TYPE = hidden NAME = borders VALUE = \"",$borders,"\">\n";

  $order_length = "";
  foreach $key (keys %order_length){
    $order_length .= "$key=$order_length{$key}\n";
  }
  print GRAPH "<INPUT TYPE = hidden NAME = order_length VALUE = \"",$order_length,"\">\n";

  print GRAPH "<P>You can redraw the graph with new cut-off score: ";
  print GRAPH "<INPUT TYPE = text NAME = SCORE VALUE = $properties{SCORE} SIZE = 3> ";
  print GRAPH "<INPUT TYPE = submit NAME = submit VALUE = Redraw><BR>\n";

  print GRAPH "<P>Areas underlined with red on this graph represent conserved regions. They are <B>clickable</B> and connected with lists of overlapping clusters of conserved motifs and their sequences.<BR>\n";
  print GRAPH "<IMG SRC = $graph_link ALT = 'Please wait for this image to load.' USEMAP = #Map BORDER = 0><BR>\n";
  print GRAPH "<MAP NAME = Map>\n";
  foreach (@enh){
    /(\d+)-(\d+)/;
    $start = $1 ; $end = $2;
    $coord1 = $marg+($start-$properties{'REF_SEQ_START'})/$scale;
    $coord2 = $marg+$gc_ruler+$conservation_ruler;
    $coord3 = $marg+($end-$properties{'REF_SEQ_START'})/$scale;
    $coord4 = $marg+$gc_ruler+$conservation_ruler+$cluster_ruler;
    $coord = "$coord1,$coord2,$coord3,$coord4\n";
    $sub_clusters = "$properties{URL}"."cgi-bin/sub_clusters.pl?START=$start&END=$end&FILE=$clusters_file&REF=$properties{$properties{'REF_SPECIES'}}&SPECIES=$species";
    print GRAPH "<area shape=rect coords=$coord href=$sub_clusters BORDER = 0>\n";
  }
  print GRAPH "</MAP>\n";
  print GRAPH "</FORM></BODY></HTML><BR>\n";

  %properties;
}


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
