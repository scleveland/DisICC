package Statistics;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;

sub Average
{
	my(@array) = @_;
	my($sum) = 0;
	for(my $i = 0; $i <= scalar(@array) - 1; $i++)
	{
		$sum += $array[$i];
	}
	return ($sum/scalar(@array));
}
sub SE
{
	my($numbers, $media) = @_;
	my($Sum_squares) = 0;
	my($i);
	for($i = 0; $i <= scalar(@$numbers) - 1; $i++)
	{
		$Sum_squares += (($$numbers[$i] - $media)**2);
	}
	my($var) = sqrt($Sum_squares/($i + 1));
	return $var;
}
sub Z
{
	my($mean, $var, $Numbers, $Z_val) = @_;
	my($i) = 0;
	while($i <= scalar(@$Numbers) - 1)
	{
		$$Z_val[$i] = ($$Numbers[$i] - $mean)/$var;
		$i += 1;
	}
}
sub Correlation
{
	my($sampleA, $sampleB) = @_;
	my(@A, @B, $meanA, $meanB) = (0,0,0,0);
	for(my $k = 0; $k <= scalar(@$sampleA) - 1; $k++)
	{
		$A[$k] = $$sampleA[$k];
	}
	for(my $k = 0; $k <= scalar(@$sampleB) - 1; $k++)
	{
		$B[$k] = $$sampleB[$k];
	}
	$meanA = Average(@A);
	$meanB = Average(@B);
	my($i) = 0;
	my($nominator) = 0;
	my($X) = 0;
	my($Y) = 0;
	my($correl);
	while($i <= scalar(@$sampleA) - 1)
	{
		if($$sampleB[$i])
		{
			$nominator += ($$sampleA[$i] - $meanA) * ($$sampleB[$i] - $meanB);
			$Y += (($$sampleB[$i] - $meanB)**2);
		}else
		{
			$nominator += 0;
			$Y += 0;
		}
		$X += (($$sampleA[$i] - $meanA)**2);
			
		$i += 1;
	}
	if(($X > 0) && ($Y > 0))
	{
		$correl = ($nominator/sqrt($X*$Y));
	}else
	{
		$correl = 0;
	}
	return $correl;
}
sub Z_probabilities
{
	my($Z_value, $Z_scores) = @_;
	my($score) = 0;
	$score += $Z_value;# 68
	my($prob) = 0;
	if($score <= 0)
	{
		$score *= (-1);
	}
	if($score == 0)
	{
		$prob = 0.51;
	}elsif($score >= 3.1)
	{
		$prob = 0.0001;
	}elsif(($score > 0) && ($score <= 3.08))
	{
		$score *= 100;
		my($score_red) = int($score);
		$prob += (1 - (($$Z_scores[$score_red] + $$Z_scores[$score_red + 1])/2));
	}
	my($prob_red) = sprintf "%4f", $prob;
	return $prob_red;
}
1;