#!/usr/bin/perl

sub trim($)
{
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

while(<STDIN>) {
  if(/UW/) {
    @head_split = split(/[\s\t]/);    
    open(PWMFILE, ">$head_split[0].pwm");
    print PWMFILE "#$head_split[0]\n";
    print PWMFILE "A	C	G	T\n";
  } elsif(/^\s*$/) {
    close(PWMFILE);
  } else {
    print PWMFILE $_;
  }
}

