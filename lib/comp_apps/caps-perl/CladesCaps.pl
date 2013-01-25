#!/usr/bin/perl -w
# This program splits inputfile into subfile with different user pre-specified clades, it runs CAPS on each file, stores the coevolving sites, compares the list of sites between all the analyses and picks comon sites among them as real coevolution.
#We provide the program by a file with the original set of sequences in fasta and with a file containing the clades that have to be taken off every run time
#the clades file has to contain the name of the sequences performing each clade in one line separated by a space

use strict;
use warnings;
use File_reader;
use Automatic_caps_module;

my(@clades, $input, $output);
$input = "Clades";
unless(open(CLADES, $input))
{
	print "Error: Could not find the file $input !!\n";
	exit;
}
###STORING THE PRE-SPECIFIED CLADES IN AN ARRAY###
@clades = <CLADES>;
my($i) = 0;
foreach my $line(@clades)
{
	chomp $line;
	$clades[$i] = $line;
	$i += 1;
}
###READING THE OPTIONS OF THE CONTROL FILE FOR CAPS###
my(@ctl_options, $control);
$control = "CAPS.ctl";
File_reader::Control_reader($control, \@ctl_options);
###RUNNING CAPS WITH ALL THE SEQUENCES INCLUDED AND STORING GROUPS OF CO-EVOLUTION IN AN ARRAY###
system './CAPS.pl';
my($file_group, @groups_of_coevolution, $outcaps, @information_out); #sites of the same group will be separated by space, each group is going to be separated by tabulation, whereas groups belonging to different running processes will be separated as elements of the array @groups_of_coevolution
$outcaps = $ctl_options[2];
open(OUT, $outcaps);
@information_out = <OUT>;
close OUT;
Automatic_caps_module::Parse_groups_coevolution(\@information_out, \$file_group);
my($pointer_to_final_groups) = 0;
$groups_of_coevolution[$pointer_to_final_groups] = $file_group;
$pointer_to_final_groups += 1;

###BUILDING FILES WITH THE SEQUENCES SUBSETS###
my(@sequences, @seq_names, @nuc_seq);
File_reader::Format_reader($ctl_options[0], \@sequences, \@seq_names, \@nuc_seq);
my($number_files) = 0;
my(@files, @array);
while($number_files <= scalar(@clades) - 1)
{
	my(@array);
	@array = split(/\s+/, $clades[$number_files]);
	$files[$number_files] = $number_files . ".in";
	open(IN, ">$files[$number_files]");
	$i = 0;
	while($i <= scalar(@sequences) - 1)
	{
		my($j) = 0;
   	 	my($name_detected) = 0;
		while($j <= scalar(@array) - 1)
		{
    		#	chomp $seq_names[$i];
			if($seq_names[$i] eq $array[$j])
			{
				$name_detected += 1;
			}
			$j += 1;
		}
		if($name_detected == 0)
		{
			print IN ">" . $seq_names[$i] . "\n" . $sequences[$i] . "\n";
		}
		$i += 1;
	}
	$number_files += 1;
}

###RUNNING CAPS AUTOMATICALLY FOR EACH ONE OF THE FILES GENERATED AND STORING GROUPS OF COEVOLUTION IDENTIFIED BY CAPS###
$i = 0;
my($out_subset);
while($i <= scalar(@files) - 1)
{
	@array = split(/in/, $files[$i]);
	$out_subset = $array[0] . "out";
	Automatic_caps_module::modify_control($files[$i], $out_subset);
	system './CAPS.pl';
	File_reader::Control_reader($control, \@ctl_options);
	$outcaps = $ctl_options[2];
	open(OUT2, $outcaps);
	@information_out = <OUT2>;
	close OUT2;
	Automatic_caps_module::Parse_groups_coevolution(\@information_out, \$file_group);
	$groups_of_coevolution[$pointer_to_final_groups] = $file_group;
	$pointer_to_final_groups += 1;
	$i += 1;
}

###FINDING THE DIFFERENT FINAL GROUPS OF COEVOLUTION###
$i = 0;
my(@array_initial_group) = split(/\t+/, $groups_of_coevolution[0]);
my($FuncStruct) = "FuncStrucPairs.csv";
open(FUNC, ">$FuncStruct");
print FUNC "\t\tDetecting functional/structural important pairs of co-evolving amino acids\n\n";
print FUNC "Group,Site1,Site2\n";

my($siteA, $siteB, $detectedA, $detectedB, $groups_initial) = (0,0,0,0,0);
while($groups_initial <= scalar(@array_initial_group) - 1)
{
	$i = 0; my($j) = 1;
	my(@array_for_this_group) = split(/\s+/, $array_initial_group[$groups_initial]);
	while($i <= scalar(@array_for_this_group) - 2)
	{
		while($j <= scalar(@array_for_this_group) - 1)
		{
			$detectedA = 0; $detectedB = 0;
			$siteA = $array_for_this_group[$i];
			$siteB = $array_for_this_group[$j];
			my($pointer_to_filegroups) = 1;
			my($both_detected) = 0;
			while($pointer_to_filegroups <= scalar(@groups_of_coevolution) - 1)
			{
				my(@array_next_group) = split(/\t+/, $groups_of_coevolution[$pointer_to_filegroups]);
				my($screen) = 0;
				while($screen <= scalar(@array_next_group) - 1)
				{
					my(@array_for_eachgroup) = split(/\s+/, $array_next_group[$screen]);
					my($screenb) = 0;
					
					while($screenb <= scalar(@array_for_eachgroup) - 1)
					{
						if($siteA eq $array_for_eachgroup[$screenb])
						{
							$detectedA += 1;
							last;
						}
						$screenb += 1;
					}
					$screenb = 0;
					while($screenb <= scalar(@array_for_eachgroup) - 1)
					{
						if($siteB eq $array_for_eachgroup[$screenb])
						{
							$detectedB += 1;
							last;
						}
						$screenb += 1;
					}
					if(($detectedA > 0) && ($detectedB > 0))
					{
						last;
					}
					$screen += 1;
				}
				if(($detectedA > 0) && ($detectedB > 0))
				{
					$both_detected += 1;
				}
				$pointer_to_filegroups += 1;
			}
			if($both_detected == scalar(@groups_of_coevolution) - 1)
			{
				print FUNC ($groups_initial + 1) . "," . $siteA . "," . $siteB . "\n";
			}
			$j += 1;
		}
		$i += 1;
		$j = $i + 1;
	}
	$groups_initial += 1;
}

exit;
