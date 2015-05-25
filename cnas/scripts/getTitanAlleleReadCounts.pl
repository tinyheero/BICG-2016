#!/usr/bin/env perl
use strict;

my $filename = $ARGV[0];
open FILE, $filename or die $1;

while (<FILE>) {
	if ($. == 1){
		print "chr\tposition\tref\trefCount\tNref\tNrefCount\n"
	} else {
		my $line = $_;
		my @elements = split("\t", $line);
		my @values = split(",", $elements[4]);
		my $refCount = $values[0] + $values[1];
		my $non_refCount = $values[2] + $values[3];
		if ($elements[3] eq "" ){ $elements[3] = "N" }
		print "$elements[0]\t$elements[1]\t$elements[2]\t$refCount\t$elements[3]\t$non_refCount" . "\n";
	}
}
