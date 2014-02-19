#!/usr/bin/perl

sub trim($)
{
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}


while(<STDIN>) {
  @names = split(/[\s\t]/);    
  $tfname = $names[2];
  $filename = $names[0];
  open(PWMFILE, "<", "$filename.pfm");
  open(PWMOUT, ">", "../$tfname.$filename.pwm");

  @A = split(/[\s\t]+/, trim(readline(PWMFILE)));
  @C = split(/[\s\t]+/, trim(readline(PWMFILE)));
  @G = split(/[\s\t]+/, trim(readline(PWMFILE)));
  @T = split(/[\s\t]+/, trim(readline(PWMFILE)));

  my $nbases = scalar @A;
  print PWMOUT "#$tfname $nbases\n";
  print PWMOUT "A\tC\tG\tT\n";
  for($i = 0; $i < $nbases; $i++) {
    print PWMOUT "$A[$i]\t$C[$i]\t$G[$i]\t$T[$i]\n";
  }
  close(PWMOUT);
  close(PWMFILE);
}

