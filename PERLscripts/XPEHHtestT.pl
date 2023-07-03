#! /usr/bin/perl
use strict;
use warnings;
use Statistics::TTest;
my $input = shift;
my $inf = shift;
my $sub = shift;
my @all;
my @fusionregion;
open I,"<$input" or die "Usage: perl $0 input inf sub";
print STDERR "File reading............\n";
while (<I>)
{
        my ($pos,$gpos,$p1,$ihh1,$p2,$ihh2,$xpehh,$normxpehh,$crit) = split(/\t/);
        if ($normxpehh eq "NA" ) {}
        elsif ( $normxpehh >  0 )
        {
        if( $pos >= $inf - 1000000 && $gpos <= $sub + 1000000)
        {
                push @all,$normxpehh;
                push @fusionregion,$normxpehh;
        }
        else
        {
                push @all,$normxpehh;
        }
        }
        else {}
}

print STDERR "All site loaded............\n";
my $ttest = new Statistics::TTest;
#$ttest->set_significance(99);
$ttest->load_data(\@all,\@fusionregion);
$ttest->output_t_test();
print STDERR "T-test done.................\nplot start ................\n";
