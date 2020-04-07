#!/usr/bin/perl
# BamSubset.pl separates bamfiles via samtools. Faster alternative to subset-bam for a single core, and if clusters > cores. Work of Thomas Burkard 04/2020

use strict;
use Getopt::Long;

my $bam="";
my $clusterFile="";
my $clusterCol="1"; #0-based, -1 = use all BC
my $bcCol="0"; #0-based
my $sep=",";
my $tag="CB";

GetOptions('bam:s'=>\$bam,
           'cluster:s'=>\$clusterFile,
    );

open(CLU, $clusterFile);
my $bc2cl = {};
my $clCnt = {};

$clCnt->{"none"} = 0;
while (<CLU>) {
    chomp;
    my @col = split $sep;

    $clCnt->{$col[$clusterCol]}++;
    $bc2cl->{$col[$bcCol]} = $col[$clusterCol];
}
close CLU;

my $clFh = {};
foreach my $cl (keys %$clCnt) {
    open($clFh->{$cl}, "| samtools view -bS - > $cl.bam");
}

open(BAM, "samtools view -h $bam |");

while (my $line=<BAM>) {

    if ($line=~/^@/) {
        foreach my $cl (keys %$clCnt) {
            print {$clFh->{$cl}} $line;
        }
    } else {
        if (($line=~/$tag:Z:(\S+)/) && exists($bc2cl->{$1}))
        {
            print {$clFh->{$bc2cl->{$1}}} $line;
        } else {
            print {$clFh->{"none"}} $line;
        }
    }
}
close BAM;
