package CAPS_module;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;
use File_reader;
use Statistics;

sub Read_blosum
{
	my($file, $B_array, $B_aa) = @_;
	my($j, $i) = (0, 1);
	my(@array, @B_aa, $in);
	$in = File_reader::Input_file($$file);
	while(<$in>)
	{
		if($_ =~ /^(BLOSUM62)/)
		{
			next;
		}
		$_ =~ s/ {2}/ /g;
		@array = split(/\s/, $_);
		$$B_aa[$j] = $array[0];
		$$B_aa[$j] =~ s/\s+//g;
		for($i = 1; $i <= 21; $i++)
		{
			$$B_array[$j][$i -1] = $array[$i];
		}
		$j += 1;  
	}
}
sub Estimate_D
{
	my($seqs, $D, $D_corr, $rel_dist, $correction, $Blosums, $B_aa, $aa_differences) = @_;
	my($longitud, $number_Ds,) = (0,0);
	my($aa1, $aa2);
	while($longitud <= length($$seqs[0]) - 2)
	{
		my(@theta);
		my($loop1, $loop2, $number_thetas) = (0,1,0);
		if($$aa_differences[$longitud] > 0)
		{
		while($loop1 <= scalar(@$seqs) - 2)
		{
			my($i) = 0;
			$aa1 = substr($$seqs[$loop1], $longitud, 1);
			while(($aa1 ne $$B_aa[$i]) && ($i < scalar(@$B_aa) - 1))
			{
				$i += 1;
			}
			while($loop2 <= scalar(@$seqs) - 1) #Line 50
			{
				my($j) = 0;
				$aa2 = substr($$seqs[$loop2], $longitud, 1);
				while(($aa2 ne $$B_aa[$j]) && ($j < scalar(@$B_aa) - 1))
				{
					$j += 1;
				}
				$theta[$number_thetas] = $$Blosums[$i][$j];
				$number_thetas += 1;
				$loop2 += 1;
			}
			$loop1 += 1;
			$loop2 = $loop1 + 1;
		}
		my($Theta_aver) = Statistics::Average(@theta);
		for(my $k = 0; $k <= scalar(@theta) - 1; $k++)
		{
			$$D[$longitud][$k] = ($theta[$k] - $Theta_aver)**2;#square would take variances
		}
		if($$correction == 1) #Linea number 70
		{
			for(my $k = 0; $k <= scalar(@theta) - 1; $k++)
			{
				if($$rel_dist[$k]  > 0)
				{
					$$D_corr[$longitud][$k] = $$D[$longitud][$k]/$$rel_dist[$k];
				}else
				{
					$$D_corr[$longitud][$k] = $$D[$longitud][$k];
				}
			}
		}
		}else
		{
			my($number_comp) = ((scalar(@$seqs)) * (scalar(@$seqs) - 1))/2;
			for(my $k = 0; $k <= $number_comp - 1; $k++)
			{
				$$D[$longitud][$k] = 0;
				$$D_corr[$longitud][$k] = 0;
			}
		}
		$longitud += 1;
	}
}
sub Co_evolution_intra
{
	my($D, $D_correct, $Correl, $information13, $seq_number, $length_seq) = @_;
	my($i, $j, $m, $n, $number_correlations, $num) = (0,1,0,0,0,0);
	my($total_number_comp) = ((scalar(@$D) - 1) * scalar(@$D))/2;
	print "\nComputing pairwise amino acid site comparisons: \n";
	while($i <= scalar(@$D) - 2)
	{
		my(@D_column);
		$m = 0;
		while($m <= (($$seq_number * ($$seq_number - 1))/2) - 1)
		{
			if($$information13 == 0)
			{
				$D_column[$m] = $$D[$i][$m];
			}else
			{
				$D_column[$m] = $$D_correct[$i][$m];
			}
			$m += 1;
		}
		my(@D_column2);
		while($j <= scalar(@$D) - 1)
		{
			$n = 0;
			while($n <= (($$seq_number * ($$seq_number - 1))/2) - 1)
			{
				if($$information13 == 0)
				{
					$D_column2[$n] = $$D[$j][$n];
				}else
				{
					$D_column2[$n] = $$D_correct[$j][$n];
				}
				$n += 1;
			}
			$$Correl[$number_correlations] = Statistics::Correlation(\@D_column, \@D_column2);
			$number_correlations += 1;
			$j += 1;
			$num += 1;
			my($percent) = ($num/$total_number_comp) * 100;
			my($percent_red) = sprintf "%.0f", $percent;
			print "\r" . ($percent_red) . "%";
		}
		$i += 1;
		$j = $i + 1;
	}
	print "\r100% of comparisons computed\n";
}
sub Co_evolution_inter
{
	my($D, $D2, $D_correct, $D2_correct, $Correl, $information13, $seq_number) = @_;
	my($i, $j, $m, $n, $number_correlations) = (0,0,0,0,0);
	while($i <= scalar(@$D) - 1)
	{
		my(@D_column);
		$m = 0;
		while($m <= (($$seq_number * ($$seq_number - 1))/2) - 1)
		{
			if($$information13 == 0)
			{
				$D_column[$m] = $$D[$i][$m];
			}else
			{
				$D_column[$m] = $$D_correct[$i][$m];
			}
			$m += 1;
		}
		my(@D_column2);
		while($j <= scalar(@$D2) - 1)
		{
			$n = 0;
			while($n <= (($$seq_number * ($$seq_number - 1))/2) - 1)
			{
				if($$information13 == 0)
				{
					$D_column2[$n] = $$D2[$j][$n];
				}else
				{
					$D_column2[$n] = $$D2_correct[$j][$n];
				}
				$n += 1;
			}
			$$Correl[$number_correlations] = Statistics::Correlation(\@D_column, \@D_column2);
			$number_correlations += 1;
			$j += 1;
		}
		$i += 1;
		$j = 0;
	}
}
sub D_average
{
	my($D, $D_aver, $number_seq) = @_;
	my($i, $j) = (0,0);
	while($i <= scalar(@$D) - 1)
	{
		$j = 0;
		my($D_sum) = 0;
		while($j <= (($$number_seq * ($$number_seq - 1))/2) - 1)
		{
			$D_sum += $$D[$i][$j];
			$j += 1;
		}
		$$D_aver[$i] = $D_sum/$j;
		$i += 1;
	}
}
sub Get_threshold
{
	my($alpha, $pseudosamples, $value, $correlations) = @_;
	my(@numbers);
	my($i, $Mean, $SE) = (0,0,0);
	my($random_number);
	while($i < $$pseudosamples)
	{
		$numbers[$i] = $$correlations[int(rand(scalar(@$correlations)))];
		$i += 1;
	}
	$Mean = Statistics::Average(@numbers);
	$SE = Statistics::SE(\@numbers, $Mean);
	if($$alpha == 0.05)
	{
		 $$value = ((1.69 * $SE) + $Mean);
	}
	if(($$alpha < 0.05) && ($$alpha > 0.025))
	{
		$$value = ((1.95 * $SE) + $Mean);
	}
	if(($$alpha <= 0.025) && ($$alpha > 0.01))
	{
		$$value = ((1.97 * $SE) + $Mean);
	}
	if(($$alpha <= 0.01) && ($$alpha > 0.001))
	{
		$$value = ((2.33 * $SE) + $Mean);
	}
	if($$alpha <= 0.001)
	{
		$$value = ((3.1 * $SE) + $Mean);
	}
}			
1;