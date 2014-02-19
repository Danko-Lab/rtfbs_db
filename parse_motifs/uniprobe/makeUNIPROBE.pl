#!/usr/bin/perl

sub trim($)
{
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}


while(<STDIN>) {
  @names = split(/[\s\t+]/);
  my $lenNames = scalar @names;

  ## Get the TF name.
  if( $lenNames > 1 ) {
   $tfname = $names[1];
   $iupac  = $names[5];
  } else {
   $tfname = $names[0];
   <STDIN>;
  }

  ## Get ATCG.
  @A = split(/[\s\t]+/, trim(<STDIN>));
  @C = split(/[\s\t]+/, trim(<STDIN>));
  @G = split(/[\s\t]+/, trim(<STDIN>));
  @T = split(/[\s\t]+/, trim(<STDIN>));

  my $nbases = scalar @A;
  print "#$tfname $iupac $nbases\n"; #PWMOUT
  print "A\tC\tG\tT\n";
  for($i = 1; $i < $nbases; $i++) {
    print "$A[$i]\t$C[$i]\t$G[$i]\t$T[$i]\n";
  }
#  close(PWMOUT);
}

