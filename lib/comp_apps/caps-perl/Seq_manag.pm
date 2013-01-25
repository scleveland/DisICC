package Seq_manag;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;
my(%codigo_genetico) = (
	'TTT' => 'F', 
	'TTC' => 'F',
	'TTA' => 'L',
	'TTG' => 'L',
	'CTT' => 'L',
	'CTC' => 'L',
	'CTA' => 'L',
	'CTG' => 'L',
	'ATT' => 'I',
	'ATC' => 'I',
	'ATA' => 'I',
	'ATG' => 'M',
	'GTT' => 'V',
	'GTC' => 'V',
	'GTA' => 'V',
	'GTG' => 'V',
	'TCT' => 'S',
	'TCC' => 'S',
	'TCA' => 'S',
	'TCG' => 'S',
	'CCT' => 'P',
	'CCC' => 'P',
	'CCA' => 'P',
	'CCG' => 'P',
	'ACT' => 'T',
	'ACC' => 'T',
	'ACA' => 'T',
	'ACG' => 'T',
	'GCT' => 'A',
	'GCC' => 'A',
	'GCA' => 'A',
	'GCG' => 'A',
	'TAT' => 'Y',
	'TAC' => 'Y',
	'CAT' => 'H',
	'CAC' => 'H',
	'CAA' => 'Q',
	'CAG' => 'Q',
	'AAT' => 'N',
	'AAC' => 'N',
	'AAA' => 'K',
	'AAG' => 'K',
	'GAT' => 'D',
	'GAC' => 'D',
	'GAA' => 'E',
	'GAG' => 'E',
	'TGT' => 'C',
	'TGC' => 'C',
	'TGG' => 'W',
	'CGT' => 'R',
	'CGC' => 'R',
	'CGA' => 'R',
	'CGG' => 'R',
	'AGT' => 'S',
	'AGC' => 'S',
	'AGA' => 'R',
	'AGG' => 'R',
	'GGT' => 'G',
	'GGC' => 'G',
	'GGA' => 'G',
	'GGG' => 'G',
	'TAA' => '-',
	'TAG' => '-',
	'TGA' => '-',
	'---'  => '-',);
	
sub codon_aminoacido
{
	my($codon) = @_;
	$codon = uc $codon;###si esta en minúscula lo transforma a mayúscula

if(exists $codigo_genetico{$codon})
	{
		return $codigo_genetico{$codon};
	}else
	{
	print STDERR "Bad codon \"$codon\"!!\n";
	exit;
	}
}
####Subrutina para la traducción de una secuencia de aminoácidos a una de proteínas###
sub nuc_amino
{
	my($dna) = @_;
	my($protein) = '';
	for(my $i = 0; $i <= (length($dna) - 2); $i += 3)
	{
		$protein .= codon_aminoacido(substr($dna, $i, 3));
	}
	return $protein;
}
1;