#! /usr/bin/perl
use strict;
use warnings;
my $dir = shift;
my $output = shift;
open O, ">$output" or die "perl $0 dir output";
print O "chr\ttargetcount\trefcount\n";
my @input = glob("$dir/*.xpehh.out.norm");
for my $input(@input)
{
open I, "<$input" or die "perl $0 dir output";
my $targetcount = 0;
my $refcount = 0;
while(<I>)
{
		my ($id,$pos,$gpos,$p1,$ihh1,$p2,$ihh2,$xpehh,$normxpehh,$crit) = split(/\t/) ;
		if ( $normxpehh >= 2 )
		{
			$targetcount = $targetcount +1;
		}
		elsif ( $normxpehh <= -2 )
		{
			$refcount = $refcount + 1 ;
		}
		else{}
}
print O "$input\t$targetcount\t$refcount\n";
$targetcount = 0;
$refcount = 0;
}
