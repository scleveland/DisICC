package CompMut;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;
use Statistics;

###Subroutine to calculate the hydrophobicity for each site at each coevolution group###
sub Hydrophobicity
{
	my($sequences, $coevolvings, $hydrophobicities, $hydros_for_site, $MW) = @_;
	my($i, $j) = (0,0);
	my(@coevols);
	while($i <= scalar(@$coevolvings) - 1)
	{
		my(@coevol);
		@coevol = split(/\s+/, $$coevolvings[$i]);
		$j = 0;
		while($j < scalar(@coevol))
		{
			my(@coevol_split);
			@coevol_split = split(/\(/, $coevol[$j]);
			if($j == 0)
			{
				$coevols[$i] = $coevol_split[0];
			}else
			{
				$coevols[$i] .= " " . $coevol_split[0];
			}
			$j += 1;
		}
		$i += 1;
	}
	$i = 0; $j = 0;
	my($bucle1, $bucle2) = (0,0);
	while($bucle1 <= scalar(@coevols) - 1)
	{
		$bucle2 = 0;
		my(@array) = split(/\s+/, $coevols[$bucle1]);
		while($bucle2 <= scalar(@array) - 1)
		{
			$i = 0;
			$$hydros_for_site[$bucle1][$bucle2] = '';
			$$MW[$bucle1][$bucle2] = '';
			while($i <= scalar(@$sequences) - 1)
			{
				my($aa);
				my($pos) = $array[$bucle2] - 1;
				$aa = substr($$sequences[$i], $pos, 1);
				my($search) = 0;
				while($search <= scalar(@$hydrophobicities) - 1)
				{
					my(@array2);
					@array2 = split(/\t+/, $$hydrophobicities[$search]);
					if($aa eq $array2[1])
					{
						if(length($$hydros_for_site[$bucle1][$bucle2]) > 0)
						{
							$$hydros_for_site[$bucle1][$bucle2] .= ' '  . $array2[2];
							$$MW[$bucle1][$bucle2] .= ' ' . $array2[4];
						}else
						{
							$$hydros_for_site[$bucle1][$bucle2] = $array2[2];
							$$MW[$bucle1][$bucle2] = $array2[4];
						}
						last;
					}
					if($aa eq "-")
					{
						if(length($$hydros_for_site[$bucle1][$bucle2]) > 0)
						{
							$$hydros_for_site[$bucle1][$bucle2] .= ' '  . 0;
							$$MW[$bucle1][$bucle2] .= ' ' . 0;
						}else
						{
							$$hydros_for_site[$bucle1][$bucle2] = 0;
							$$MW[$bucle1][$bucle2] = 0;
						}
						last;
					}
					$search += 1;
				}
				$i += 1;
			}
			$bucle2 += 1;
			
		}
		$bucle1 += 1;	
	}
}

###Subroutine to calculate the variability in the hydrophobicity at each site within each coevolution group###
sub HydroVar
{
	my($Hydros_per_site, $num_group_coevol, $num_seqs, $coevol, $Hv, $MW, $MWv) = @_;
	my($i) = 0;
	while($i <= $$num_group_coevol - 1)
	{
		my(@array_each_group) = split(/\s+/, $$coevol[$i]);
		my($j) = 0;
		while($j <= scalar(@array_each_group) - 1)
		{
			my($bucle1, $bucle2, $k) = (0,1,0);
			$$Hv[$i][$j] = '';
			$$MWv[$i][$j] = '';
			my(@array_each_site) = split(/\s+/, $$Hydros_per_site[$i][$j]);
			my(@array_site_MW) = split(/\s+/, $$MW[$i][$j]);
			while($bucle1 <= $$num_seqs - 2)
			{
				while($bucle2 <= $$num_seqs - 1)
				{
					if(length($$Hv[$i][$j]) == 0)
					{
						$$Hv[$i][$j] = sqrt(($array_each_site[$bucle1] - $array_each_site[$bucle2])**2);
						$$MWv[$i][$j] = sqrt(($array_site_MW[$bucle1] - $array_site_MW[$bucle2])**2);
					}else
					{
						if(($array_each_site[$bucle1]) && ($array_each_site[$bucle2]))
						{
							$$Hv[$i][$j] .= " " . sqrt(($array_each_site[$bucle1] - $array_each_site[$bucle2])**2);
							$$MWv[$i][$j] .= " " . sqrt(($array_site_MW[$bucle1] - $array_site_MW[$bucle2])**2);
						}else
						{
							$$Hv[$i][$j] .= " " . 0;
							$$MWv[$i][$j] .= " " . 0;
						}
					}
					$bucle2 += 1;
				}
				$bucle1 += 1;
				$bucle2 = $bucle1 + 1;
			}
			$j += 1;
		}
		$i += 1;
	}
}

##Subroutine to estimate the correlation in hydrophobicity between two sites belonguing to the same group of coevolution###
sub HydroCorr
{
	my($groups_of_coevolution, $HydCorr, $Hv, $MWCorr, $MWv) = @_;
	my($i) = 0;
	while($i <= scalar(@$groups_of_coevolution) - 1)
	{
		my(@array) = split(/\s+/, $$groups_of_coevolution[$i]);
		my($bucle1, $bucle2, $par) = (0,1,0);
		while($bucle1 <= scalar(@array) - 2)
		{
			while($bucle2 <= scalar(@array) - 1)
			{
				my(@array_site1) = split(/\s+/, $$Hv[$i][$bucle1]);
				my(@array_site2) = split(/\s+/, $$Hv[$i][$bucle2]);
				my(@arrayMW_site1) = split(/\s+/, $$MWv[$i][$bucle1]);
				my(@arrayMW_site2) = split(/\s+/, $$MWv[$i][$bucle2]);
				$$HydCorr[$i][$par] = Statistics::Correlation(\@array_site1, \@array_site2);
				$$MWCorr[$i][$par] = Statistics::Correlation(\@arrayMW_site1, \@arrayMW_site2);
				$par += 1;
				$bucle2 += 1;
			}
			$bucle1 += 1;
			$bucle2 = $bucle1 + 1;
		}
		$i += 1;
	}
}
1;