package PerlModule;
require Exporter;
@ISA = (Exporter);

@Exporter = qw (out);

use strict;
use warnings;
#subroutine that takes as inputs a sequence and a number indicating number of mutations and returns the mutated sequence
sub mutate
{
	my($sequence, $number) = @_;
	my(@nucleotides) = ('A','C','G','T');
	my $j = 0;
	while($j <= $number)
	{
		my $random_position = random_pos($sequence);
		my $random_nucleotide = random_nuc(@nucleotides);
		substr($sequence,$random_position,1,$random_nucleotide);
		$j += 1;
	}
	return $sequence;
}
# Subroutine that returns the random position between 0 and the length of the passed sequence
sub random_pos
{
	my($sequence) = @_;
	my($position) = 0;
	$position = int(rand(length($sequence)));
	return $position;
}
# Subroutine to return a random nucleotide
sub random_nuc
{
	my(@nucleotides) = @_;
	my($nuc);
	$nuc = $nucleotides[int(rand(scalar @nucleotides))];
	return $nuc;
}
#subroutine that selects a random sequence among a set of provided sequences
sub random_seqs
{
	my(@seqs) = @_;
	my($selected_sequence) = $seqs[int(rand(scalar @seqs))];
	return $selected_sequence;
}
#subroutine that calculates the number of codons from each codon type
sub Codon_usage
{
	my($sequence, $codon_frequencies) = @_;
	my $i = 0;
	my(@codon_table) = ('TTT',
	'TTC',
	'TTA',
	'TTG',
	'CTT',
	'CTC',
	'CTA',
	'CTG',
	'ATT',
	'ATC',
	'ATA',
	'ATG',
	'GTT',
	'GTC',
	'GTA',
	'GTG',
	'TCT',
	'TCC',
	'TCA',
	'TCG',
	'CCT',
	'CCC',
	'CCA',
	'CCG',
	'ACT',
	'ACC',
	'ACA',
	'ACG',
	'GCT',
	'GCC',
	'GCA',
	'GCG',
	'TAT',
	'TAC',
	'CAT',
	'CAC',
	'CAA',
	'CAG',
	'AAT',
	'AAC',
	'AAA',
	'AAG',
	'GAT',
	'GAC',
	'GAA',
	'GAG',
	'TGT',
	'TGC',
	'TGG',
	'CGT',
	'CGC',
	'CGA',
	'CGG',
	'AGT',
	'AGC',
	'AGA',
	'AGG',
	'GGT',
	'GGC',
	'GGA',
	'GGG',
	'TAA',
	'TAG',
	'TGA',
	);
	while(my $codon = substr($sequence, $i, 3))
	{
		for(my $j = 0; $j <= 63; $j++)
		{
			if($codon eq $codon_table[$j])
			{
				$$codon_frequencies[$j] += 1;
			}else
			{
				$$codon_frequencies[$j] += 0;
			}
		}
		$i += 3;
	}
	return @$codon_frequencies;
}


1;
