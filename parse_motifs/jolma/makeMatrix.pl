#!/usr/bin/perl

sub trim($)
{
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

open (TFDBFILE, ">@ARGV[0]/tf.tsv");

while(<STDIN>) {
  $header = trim($_);
  print TFDBFILE "$header\n";

  @head_split = split(/[\s\t]/, $header);
  open(PWMFILE, ">@ARGV[0]/$head_split[0].$head_split[5].pwm");

  $line = trim(<STDIN>);
  @A = split(/[\s\t]/, $line);
  $line = trim(<STDIN>);
  @C = split(/[\s\t]/, $line);
  $line = trim(<STDIN>);
  @G = split(/[\s\t]/, $line);
  $line = trim(<STDIN>);
  @T = split(/[\s\t]/, $line);

  $n = @A;
  print PWMFILE "#$header\n";
  for($i=0; $i<$n; $i++) {
    print PWMFILE "$A[$i]\t$C[$i]\t$G[$i]\t$T[$i]\n";
  }
  close(PWMFILE);
}
close(TFDBFILE);
