package Distances;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;
my %cero_fold = (
	'ATG' => 'M',
	'TGG' => 'W',
	'TAA' => '-',
	'TAG' => '-',
	'TGA' => '-'
	);
my %two_fold = (
	'TTT' => 'F',
	'TTC' => 'F',
	'TAT' => 'Y',
	'TAC' => 'Y',
	'CAT' => 'H',
	'CAC' => 'H',
	'CAA' => 'Q',
	'CAG' => 'Q',
	'AAA' => 'K',
	'AAG' => 'K',
	'GAT' => 'D',
	'GAC' => 'D',
	'GAA' => 'E',
	'GAG' => 'E',
	'AAT' => 'N',
	'AAC' => 'N',
	'TGT' => 'C',
	'TGC' => 'C');
my %four_fold = (
	'ATT' => 'I',
	'ATC' => 'I',
	'ATA' => 'I',
	'GTA' => 'V',
	'GTT' => 'V',
	'GTC' => 'V',
	'GTG' => 'V',
	'CCA' => 'P',
	'CCT' => 'P',
	'CCC' => 'P',
	'CCG' => 'P',
	'ACA' => 'T',
	'ACT' => 'T',
	'ACC' => 'T',
	'ACG' => 'T',
	'GCA' => 'A',
	'GCT' => 'A',
	'GCC' => 'A',
	'GCG' => 'A',
	'GGA' => 'G',
	'GGT' => 'G',
	'GGC' => 'G',
	'GGG' => 'G');
my %six_fold = (
	'TTA' => 'L',
	'TTG' => 'L',
	'CTA' => 'L',
	'CTT' => 'L',
	'CTC' => 'L',
	'CTG' => 'L',
	'TCA' => 'S',
	'TCT' => 'S',
	'TCC' => 'S',
	'TCG' => 'S',
	'AGT' => 'S',
	'AGC' => 'S',
	'AGA' => 'R',
	'AGG' => 'R',
	'CGA' => 'R',
	'CGT' => 'R',
	'CGC' => 'R',
	'CGG' => 'R');
# Li's method to estimate synonymous distances
my %purines = (
	'A' => 'Pur1',
	'G' => 'Pur2');
my %pyrimidines = (
	'T' => 'Pyr1',
	'C' => 'Pyr2');
sub Li_synonymous
{
	my($seq1, $seq2) = @_;
	my($nucleot) = 0;
	my(@L);
	my(@Ts);
	my(@Tv);
	my(@A);
	my(@B);
	my($codon1);
	my($codon2);
	for(my $n = 0; $n <= 2; $n++)
	{
		$L[$n] = 0;
		$Ts[$n] = 0;
		$Tv[$n] = 0;
		$A[$n] = 0;
		$B[$n] = 0;
	}
	my($longi) = 0;
	while($longi <= length($seq1) - 1)
	{
		$codon1 = substr($seq1, $longi, 3);
		$codon2 = substr($seq2, $longi, 3);
		deg_transvers($codon1, $codon2, \@L, \@Ts, \@Tv);
		$longi += 3;
	}
	Calcula_A_B(\@Ts, \@Tv, \@L, \@A, \@B);
	my($ds);
	if(($L[1] + $L[2]) > 0)
	{
		$ds += (((($L[1] * $A[1]) + ($L[2] * $A[2]))/($L[1] + $L[2])) + $B[2]);
	}else
	{
		$ds += $B[2];
	}
	return $ds;
}
sub deg_transvers
{
	my($codon1, $codon2, $L, $Ts, $Tv) = @_;
	my($nucA);
	my($nucB);
	if(defined $cero_fold{$codon1})
	{
		$$L[0] += 3;
		for(my $j = 0; $j <= 2; $j++)
		{
			$nucA = substr($codon1, $j, 1);
			$nucB = substr($codon2, $j, 1);
			if($nucA ne $nucB)
			{
				if(((defined $purines{$nucA}) && (defined $purines{$nucB})) || ((defined $pyrimidines{$nucA}) && (defined $pyrimidines{$nucB})))
				{
					$$Ts[0] += 1;
				}elsif(((defined $purines{$nucA}) && (defined $pyrimidines{$nucB})) || ((defined $pyrimidines{$nucA}) && (defined $purines{$nucB})))
				{
					$$Tv[0] += 1;
				}
			}
		}
	}elsif(defined $two_fold{$codon1})
	{
		$$L[0] += 2;
		$$L[1] += 1;
		for(my $j = 0; $j <= 2; $j++)
		{
			$nucA = substr($codon1, $j, 1);
			$nucB = substr($codon2, $j, 1);
			if($nucA ne $nucB)
			{
				if(((defined $purines{$nucA}) && (defined $purines{$nucB})) || ((defined $pyrimidines{$nucA}) && (defined $pyrimidines{$nucB})))
				{
					if($j <= 1)
					{
						$$Ts[0] += 1;
					}else
					{
						$$Ts[1] += 1;
					} 
				}elsif(((defined $purines{$nucA}) && (defined $pyrimidines{$nucB})) || ((defined $pyrimidines{$nucA}) && (defined $purines{$nucB})))
				{
					if($j <= 1)
					{
						$$Tv[0] += 1;
					}else
					{
						$$Tv[1] += 1;
					}
				}
			}
		}
	}elsif(defined $four_fold{$codon1})
	{
		$$L[0] += 2;
		$$L[2] += 1;
		for(my $j = 0; $j <= 2; $j++)
		{
			$nucA = substr($codon1, $j, 1);
			$nucB = substr($codon2, $j, 1);
			if($nucA ne $nucB)
			{
				if(((defined $purines{$nucA}) && (defined $purines{$nucB})) || ((defined $pyrimidines{$nucA}) && (defined $pyrimidines{$nucB})))
				{
					if($j <= 1)
					{
						$$Ts[0] += 1;
					}else
					{
						$$Ts[2] += 1;
					} 
				}elsif(((defined $purines{$nucA}) && (defined $pyrimidines{$nucB})) || ((defined $pyrimidines{$nucA}) && (defined $purines{$nucB})))
				{
					if($j <= 1)
					{
						$$Tv[0] += 1;
					}else
					{
						$$Tv[2] += 1;
					}
				}
			}
		}
	}elsif(defined $six_fold{$codon1})
	{
		$$L[0] += 1;
		$$L[1] += 1;
		$$L[2] += 1;
		for(my $j = 0; $j <= 2; $j++)
		{
			$nucA = substr($codon1, $j, 1);
			$nucB = substr($codon2, $j, 1);
			if($nucA ne $nucB)
			{
				if(((defined $purines{$nucA}) && (defined $purines{$nucB})) || ((defined $pyrimidines{$nucA}) && (defined $pyrimidines{$nucB})))
				{
					if($j == 0)
					{
						$$Ts[1] += 1;
					}elsif($j == 1)
					{
						$$Ts[0] += 1;
					}elsif($j == 2)
					{
						$$Ts[2] += 1;
					}					
				}elsif(((defined $purines{$nucA}) && (defined $pyrimidines{$nucB})) || ((defined $pyrimidines{$nucA}) && (defined $purines{$nucB})))
				{
					if($j == 0)
					{
						$$Tv[1] += 1;
					}elsif($j == 1)
					{
						$$Tv[0] += 1;
					}elsif($j == 2)
					{
						$$Tv[2] += 1;
					}		
				}
			}
		}
	}elsif($codon1 eq '---')
	{
		for(my $n = 0; $n <= 2; $n++)
		{
			$$L[$n] += 0;
			$$Ts[$n] += 0;
			$$Tv[$n] += 0;
		}
	}
}
				

sub Calcula_A_B
{
	my($Ts, $Tv, $L, $A, $B) = @_;
	my(@P);
	my(@Q);
	my($n) = 0;
	for($n = 0; $n <= 2; $n++)
	{
		$P[$n] = 0;
		$Q[$n] = 0;
	}
	for($n = 0; $n <= 2; $n++)
	{
		if($$L[$n] > 0)
		{
			$P[$n] += ($$Ts[$n]/$$L[$n]);
			$Q[$n] += ($$Tv[$n]/$$L[$n]);
		}else
		{
			$P[$n] += 0;
			$Q[$n] += 0;
		}
		my($razon_A) = 0;
		my($razon_B) = 0;
		$razon_A += (1 - $Q[$n] - (2 * $P[$n]));
		$razon_B += (1 - (2 * $Q[$n]));
		if($razon_A < 0)
		{
			$razon_A *= (-1);
		}
		if($razon_B < 0)
		{
			$razon_B *= (-1);
		}
		my($dividendo_A) = 0;
		my($dividendo_B) = 0;
		if($razon_A > 0)
		{
			$dividendo_A += (1/$razon_A);
		}else
		{
			$dividendo_A += 1;
		}
		if($razon_B > 0)
		{
			$dividendo_B += (1/$razon_B);
		}else
		{
			$dividendo_B += 1;
		}
		$$A[$n] += ((0.5 * log($dividendo_A)) - (0.25*log($dividendo_B)));
		$$B[$n] += (0.5 * log($dividendo_B));
	}
}
sub Average
{
	my(@numbers) = @_;
	my($sum) = 0;
	my($mean) = 0;
	for(my $i = 0; $i <= scalar(@numbers) - 1; $i++)
	{
		$sum += $numbers[$i];
	}
	$mean += ($sum/scalar(@numbers));
	return $mean;
}
#Atomic distance between two amino acids
sub atom_dist
{
	my($coordinatesA, $coordinatesB) = @_;
	my($module) = 0;
	$module += sqrt((($$coordinatesA[0] - $$coordinatesB[0])**2) + (($$coordinatesA[1] - $$coordinatesB[1])**2) + (($$coordinatesA[2] - $$coordinatesB[2])**2));
	return $module;
}
sub Poisson_dist
{
	my($seq1, $seq2) = @_;
	my($i) = 0;
	my($dist) = 0;
	my($gap) = 0;
	while($i <= length($seq1) - 1)
	{
		my($aa1) = substr($seq1, $i, 1);
		my($aa2) = substr($seq2, $i, 1); 
		if(($aa1 ne $aa2) && ($aa1 ne '-') && ($aa2 ne '-'))
		{
			$dist += 1;
		}
		if(($aa1 eq '-') || ($aa2 eq '-'))
		{
			$gap += 1;
		}
		$i += 1;
	}
	return Poisson($dist, (length($seq1) - $gap));
}		
sub Poisson
{
	my($distance, $long) = @_;
	my($poisson_d);
	if((1 - ($distance/$long)) != 0)
	{
		$poisson_d = -log(1 - ($distance/$long));
	}else
	{
		$poisson_d = 0;
	}
	return $poisson_d;
}
sub Relative_distance
{
	my($distances, $dist_sort, $Rel_dist) = @_;
	my(@relative_distance);
	for(my $i = 0; $i < scalar(@$distances); $i++)
	{
		if($$dist_sort[0] != 0)
		{
			$relative_distance[$i] = $$distances[$i]/$$dist_sort[0];
		}else
		{
			$$dist_sort[0] = 1;
			$relative_distance[$i] = $$distances[$i]/$$dist_sort[0];
		}
		$$Rel_dist[$i] = sprintf"%.4f", $relative_distance[$i];
	}
}
sub Diff_aa_column
{
	my($sequences, $array) = @_;
	my($i, $j) = (0,0);
	while($i <= length($$sequences[0]) - 1)
	{
		$j = 0;
		my(@aa_array);
		 $$array[$i] = 0;
		while($j <= scalar(@$sequences) - 1)
		{
			$aa_array[$j] = substr($$sequences[$j], $i, 1);
			if($j > 0)
			{
				my($find, $search) = (0,0);
				while($search <= scalar(@aa_array) - 2)
				{
					if($aa_array[$j] eq $aa_array[$search])
					{
						$find += 1;
						last;
					}
					$search += 1;
				}
				if($find == 0)
				{
					$$array[$i] += 1;
				}
			}
			$j += 1;
		}
		$i += 1;
	}
}		
1;