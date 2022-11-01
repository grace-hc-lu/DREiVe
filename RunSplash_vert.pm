package RunSplash_vert;
use Exporter;
use Cwd;
@ISA = ('Exporter');
@EXPORT = ('run_splash');

sub run_splash{
  my $properties = shift;
  my %properties = %$properties;
  my @min_tokens = split(/:/,$properties{'MIN_TOKENS'});
  my @min_tokens_sorted = sort {$a <=> $b} (@min_tokens);
  my $min_tokens = $min_tokens_sorted[0];
  my $density_window = $properties{'DENSITY_WINDOW'};
  my $density_tokens = $properties{'DENSITY_TOKENS'};
  my $fasta = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."fasta";  
  my $ps = "$properties{'DIR'}"."htdocs/tmp/$properties{'TITLE'}"."ps";
  my $seq_number = 0;
  if ($properties{'SELECT'} eq 'all'){
    my $grep =  `/bin/grep -i ">" $fasta`;
    my @grep = split(/\n/,$grep);
    $seq_number = scalar @grep;
  }elsif ($properties{'SELECT'} eq 'select'){
      $seq_number = $properties{'MINSPECIES'};
  }
  $SIG{ALRM} = sub{die "timeout"};
  eval{
    alarm(3600);
    $pid = fork();
    unless ($pid){  
      system "$properties{'DIR_ABS'}"."cgi-bin/Splash/splash -P standalone -a regular -q dna -i -j $seq_number -k $density_tokens -l $min_tokens -w $density_window -v -u -x 1000000 $fasta < /dev/null > /dev/null 2>&1";   
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
    $properties{'MIN_PAT_LENGTH'} = $min_tokens;
  } 
  %properties;
}  
  

