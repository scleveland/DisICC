package Automatic_caps_module;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;
use File_reader;

#Subroutine to modify the control_file of CAPS
sub modify_control
{
	my($in_data, $out_data) = @_;
	my($input, $output);
	$input = "CladesCaps.ctl";
	my($in) = File_reader::Input_file($input);
	$output = "CAPS.ctl";
	my($out) = File_reader::Output_file($output);
	my($i) = 0;
	while(<$in>)
	{
		if($i == 0)
		{
			my(@array1) = split(/\t+/, $_);
			my(@array2) = split(/:/, $array1[0]);
			print $out $array2[0] . ": " . $in_data . "\t" . $array1[1];
		}elsif($i == 1)
		{
			my(@array1) = split(/\t+/, $_);
			my(@array2) = split(/:/, $array1[0]);
			print $out $array2[0] . ": " . $in_data . "\t" . $array1[1];
		}
		elsif($i == 2)
		{
			my(@array1) = split(/\t+/, $_);
			my(@array2) = split(/:/, $array1[0]);
			print $out $array2[0] . ": " . $out_data . "\t" . $array1[1];
		}else
		{
			print $out $_;
		}
		$i += 1;
	}
}
###THIS SUBROUTINE WILL PARSE THE INFORMATION OF THE GROUPS OF COEVOLUTION GENERATED FROM CAPS
sub Parse_groups_coevolution
{
	my($information, $groups) = @_;
	my($i, $pointer, $array_pointer, $j) = (0,0,0,0);
	my(@first_array);
	###SEPARATING GROUPS OF COEVOLUTION FROM THE REST OF THE CAPS OUTPUT FILE###
	my($size) = scalar(@$information);
	while($j <= $size - 1)
	{
		if($$information[$j] =~ /^(Groups)/)
		{
			$i += 1;
		}
		if($i > 0)
		{
			$pointer += 1;
		}
		if(($pointer >= 2) && ($$information[$j] =~ /^G/))
		{
			$first_array[$array_pointer] = $$information[$j];
			$array_pointer += 1;
		}
		my($letter);
		$letter = substr($$information[$j], 0, 1);
		if(($pointer > 2) && ($letter ne 'G'))
		{
			last;
		}
		$j += 1;
	}
	###STORING GROUPS OF COEVOLUTION APPROPRIATELY
	$i = 0;
	while($i <= scalar(@first_array) - 1)
	{
		chomp $first_array[$i];
		my(@array) = split(/:/, $first_array[$i]);
		$array[1] =~ s/\s//;
		if($i == 0)
		{
			$$groups = $array[1];
		}else
		{
			$$groups .= "\t" . $array[1];
		}
		$i += 1;
	}
}
1;
