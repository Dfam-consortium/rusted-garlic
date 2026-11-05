#!/usr/bin/perl
# Per-record complexity values from SeqComplex.pm, both code paths:
#   s_<m>: per-measure subroutines (profileComplexSeq.pl / runAllMethods)
#   i_<m>: runAllMethodsInline (gatherStats.pl)
# Usage: SEQCOMPLEX_HOME=<dir with SeqComplex.pm> perl perl_complexity.pl windows.fa > perl.tsv
#        (window = whole record)
use strict; use warnings;
use File::Temp qw(tempdir);
BEGIN { die "set SEQCOMPLEX_HOME to the directory containing SeqComplex.pm\n" unless $ENV{SEQCOMPLEX_HOME} }
use lib $ENV{SEQCOMPLEX_HOME};
use SeqComplex;
my @m = qw/gc at gcs ats cpg ket pur cwf ce cz cm1 cm2 cm3 cm4 cm5 cm6 ct1 ct2 ct3 ct4 ct5 ct6 cl1 cl2 cl3 cl4 cl5 cl6/;
my (%seq, @order, $id);
open F, $ARGV[0] or die; while (<F>) { chomp; if (/^>(\S+)/) { $id=$1; push @order,$id; $seq{$id}=''; } else { $seq{$id}.=uc $_; } } close F;
chdir tempdir(CLEANUP => 1) or die;   # cz writes temp files in cwd
print join("\t","seq", map("s_$_",@m), map("i_$_",@m)), "\n";
for my $id (@order) {
  my $s = $seq{$id};
  open my $null, '>', '/dev/null'; my $old = select $null;   # ce() prints debug lines
  # Both paths die on some very short records (subroutines: log of 0 on 1 bp;
  # inline: division by zero on <= 7 bp); report NA for those.
  my %a = eval { runAllMethods($s) };
  my %b = eval { runAllMethodsInline($s) };
  select $old;
  my @row = ($id);
  for my $x (@m) { my $v=$a{$x}[0]; push @row, defined $v ? $v : 'NA'; }
  for my $x (@m) { my $v=$b{$x}[0]; push @row, defined $v ? $v : 'NA'; }
  print join("\t",@row),"\n";
}
