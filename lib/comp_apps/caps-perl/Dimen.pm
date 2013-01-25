package Dimen;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;
use File_reader;

sub parse_3D
{
	my($file, $atoms, $atom_interval1, $atom_interval2) = @_;
	my($atom);
	my($AAs) = 0;
	$atom = File_reader::Input_file($$file);
	while(<$atom>)
	{
		if($_ =~ /^(ATOM)/)
		{
			my(@array) = split(/\s+/, $_);
			if(($array[1] >= $$atom_interval1) && ($array[1] <= $$atom_interval2))
			{
				$$atoms[$AAs] = $_;
				$AAs += 1;
			}
		}else
		{
			next;
		}
	}
}##line 30
##This subroutine stores the  average atomic coordinates for each amino acid
sub parse_3D_coordinates
{
	my($atom, $coordinates) = @_;
	my($i, $atoms, $number_aas) = (0,0,0);
	my(@array);
	@array = split(/\s+/, $$atom[0]);
	$$coordinates[0] = '';
	my($aa_anterior);
	my(@sum_coordinates) = (0,0,0);
	my($size) = scalar(@array);
	if($size >= 12)
	{
		$$coordinates[0] .= $array[5] . ' ' . $array[6] .  ' ' .$array[7] .  ' ' . $array[8];
		$aa_anterior = $array[5];
		$sum_coordinates[0] += $array[6];
		$sum_coordinates[1] += $array[7];
		$sum_coordinates[2] += $array[8];
	}else
	{
		$$coordinates[0] .= $array[4] . ' ' . $array[5] .  ' ' .$array[6] .  ' ' . $array[7];
		$aa_anterior = $array[4];
		$sum_coordinates[0] += $array[5];
		$sum_coordinates[1] += $array[6];
		$sum_coordinates[2] += $array[7];
	}	
	$i += 1;
	$number_aas += 1;
	while($i <= scalar(@$atom - 1))
	{
		@array = split(/\s+/, $$atom[$i]);
		if($array[5] == $aa_anterior)
		{
			$sum_coordinates[0] += $array[6];
			$sum_coordinates[1] += $array[7];
			$sum_coordinates[2] += $array[8];
			$number_aas += 1;
		}else
		{
			$$coordinates[$atoms] = '';
			$$coordinates[$atoms] .= $aa_anterior .  ' ' . ($sum_coordinates[0]/$number_aas) . ' ' . ($sum_coordinates[1]/$number_aas) . ' ' . ($sum_coordinates[2]/$number_aas);
			$number_aas = 0;
			$atoms += 1;
			$$coordinates[$atoms] = '';
			$$coordinates[$atoms] .= $array[5] . ' ' . $array[6] . ' ' . $array[7] . ' ' . $array[8];
			@sum_coordinates = (0,0,0);
			$sum_coordinates[0] += $array[6];
			$sum_coordinates[1] += $array[7];
			$sum_coordinates[2] += $array[8];
			$number_aas += 1;
			$aa_anterior = $array[5];
		}
		$i += 1;
	}
}
			
sub Atom_distance
{
	my($coordinates1, $coordinates2) = @_;
	my($suma, $i) = (0,0);
	while($i <= 2)
	{
		if(($$coordinates1[$i]) && ($$coordinates2[$i]))
		{
			$suma += ($$coordinates1[$i] - $$coordinates2[$i])**2;
		}else
		{
			$suma += 999;
		}
		$i += 1;
	}
	return sqrt($suma);
}
1;