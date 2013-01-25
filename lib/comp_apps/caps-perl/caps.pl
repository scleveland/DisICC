#!/usr/bin/perl -w
# This program reads fasta or Mega formatted amino acid (or codon-based) alignments and determine blocks of amino acid evolution.
# Then, the different blocks are evaluated for co-evolution
use strict;
use warnings;
use File_reader;
use Distances;
use Seq_manag;
use Statistics;
use Randomisation;
use Dimen;
use caps_module;
use CompMut;
use InterSort;
my($header) =
"\n\t\t*****************************************************************
\t\tCAPS: Co-Evolution Analysis using Protein Sequences
\t\tAuthor: Mario A. Fares
\t\tCode for Inter-protein co-evolution clustering: David McNally
\t\tEvolutionary Genetics and Bioinformatics Laboratory
\t\tDepartment of Genetics
\t\tSmurfit Institute of Genetics
\t\tUniversity of Dublin, Trinity College
\t\tMathematical Model: Fares and Travers, Genetics (2006)173: 9 - 23 
\t\t********************************************************************\n";
print $header;
###INFORMATION INPUT BY THE USER######
srand(time|$$);
my($control, @information, $infile, $out, $clase, $seq_ref, $test, $R_cutoff);#Input, output files, type of data (nt or aa) and reference seq for 3D
#$control = 'CAPS.ctl';
$control = "$ARGV[0]";
File_reader::Control_reader($control, \@information);
print "\nReading Control file options....\n";
my($ctl) = 0;
while($ctl <= @information - 1)
{
	print '=> ' . $information[$ctl] . "\n";
	$ctl += 1;
}
my($coordinate_file, @atoms, @sequences, @seq_names, @seq_names2, @sequences2, @nuc_seq1, @nuc_seq2, $number_sequences);
$out = File_reader::Output_file($information[2]);
if($information[17] == 2)
{
	File_reader::Format_reader($information[0], \@sequences, \@seq_names, \@nuc_seq1);
	File_reader::Format_reader($information[1], \@sequences2, \@seq_names2, \@nuc_seq2);
	if($information[4] == 1)
	{
		File_reader::codon2aa(\@sequences);
		File_reader::codon2aa(\@sequences2);
	}
}
if(($information[17] == 0) || ($information[17] == 1))
{
	my(@Seqs, @Seqs_names, @Seqs2, @Seqs_names2);
	File_reader::Format_reader($information[0], \@sequences, \@seq_names, \@nuc_seq1);
	File_reader::Format_reader($information[1], \@sequences2, \@seq_names2, \@nuc_seq2);
	if($information[4] == 1)
	{
		File_reader::codon2aa(\@sequences);
		File_reader::codon2aa(\@sequences2);
	}
	my($seq_sp, $seq_site, $gaps_per_column) = (0,0,0);
	my(@gaps_per_site);
	while($seq_site <= length($Seqs[0]) - 1)
	{
		$seq_sp = 0;
		$gaps_per_column = 0;
		while($seq_sp <= scalar(@Seqs_names) - 1)
		{
			my($aa) = substr($Seqs[$seq_sp], $seq_site, 1);
			if($aa eq '-')
			{
				$gaps_per_column += 1;
			}
			$seq_sp += 1;
		}
		$gaps_per_site[$seq_site] = $gaps_per_column;
		$seq_site += 1;
	}
	if($information[17] == 0)
	{
		my($sp, $site) = (0,0);
		while($sp <= scalar(@Seqs_names) - 1)
		{
			$seq_names[$sp] = $Seqs_names[$sp];
			$site = 0;
			$sequences[$sp] = '';
			while($site <= length($Seqs[$sp]) - 1)
			{
				if($gaps_per_site[$site] == 0)
				{
					$sequences[$sp] .= (substr($Seqs[$sp], $site, 1));
				}
				$site += 1;
			}
			$sp += 1;
		}
	}
	$seq_site = 0;
	while($seq_site <= length($Seqs2[0]) - 1)
	{
		$seq_sp = 0;
		$gaps_per_column = 0;
		while($seq_sp <= scalar(@Seqs_names2) - 1)
		{
			my($aa) = substr($Seqs2[$seq_sp], $seq_site, 1);
			if($aa eq '-')
			{
				$gaps_per_column += 1;
			}
			$seq_sp += 1;
		}
		$gaps_per_site[$seq_site] = $gaps_per_column;
		$seq_site += 1;
	}
	if($information[17] == 0)
	{
		my($sp, $site) = (0,0);
		while($sp <= scalar(@Seqs_names2) - 1)
		{
			$seq_names2[$sp] = $Seqs_names2[$sp];
			$site = 0;
			$sequences2[$sp] = '';
			while($site <= length($Seqs2[$sp]) - 1)
			{
				if($gaps_per_site[$site] == 0)
				{
					$sequences2[$sp] .= (substr($Seqs2[$sp], $site, 1));
				}
				$site += 1;
			}
			$sp += 1;
		}
	}
}							
		
my($AAs) = 0;
my(@atom_inf) = split(/\-/, $information[10]);
my($Mean_atom_dist, $SE_atom_dist, $random_distances1, $random_distances2, @rand_coordinates1, @rand_coordinates2,@rand_vector);
my(@threeD_coordinates);

###Storing the information of the atoms in @atoms###
if($information[6] == 0)
{
	print "\nStoring Atomic information from the 3D structure....................\n";
}
if(($information[6] == 0) && ($information[3] == 0))
{
	Dimen::parse_3D(\$information[9], \@atoms, \$atom_inf[0], \$atom_inf[1]);
	Dimen::parse_3D_coordinates(\@atoms, \@threeD_coordinates);
	my($repeat,$r, $R) = (0,0,0);
	my(@first_line) = split(/\s+/, $atoms[0]);
	my($size) = scalar(@first_line);
	while($repeat <= 1000)
	{
		$r = int(rand(scalar(@atoms)));
		@first_line = split(/\s+/, $atoms[$r]);
		my(@second_line);
		$random_distances1 = $atoms[$r];
		do
		{
			$R = int(rand(scalar(@atoms)));
			@second_line = split(/\s+/, $atoms[$R]);
			
		}until(($R != $r) && ($first_line[5] != $second_line[5]));
		$random_distances2 = $atoms[$R];
		my($linea) = 0;
		my($number_atoms1, $number_atoms2) = (0,0);
		@rand_coordinates1 = (0,0,0);
		@rand_coordinates2 = (0,0,0);
		while($linea <= scalar(@atoms) - 1)
		{
			my(@array1) = split(/\s+/, $atoms[$linea]);
			my(@array2) = split(/\s+/, $random_distances1);
			my(@array3) = split(/\s+/, $random_distances2);
			my($start) = 0;
			my($aa_real, $aa_rand1, $aa_rand2, $position_sum);
			if($size >= 12)
			{
				$aa_real = $array1[5];
				$aa_rand1 = $array2[5];
				$aa_rand2 = $array3[5];
				$position_sum = 6;
			}else
			{
				$aa_real = $array1[4];
				$aa_rand1 = $array2[4];
				$aa_rand2 = $array3[4];
				$position_sum = 5;
			}
			if($aa_real == $aa_rand1)
			{
				for(my $suma = 0; $suma <= 2; $suma++)
				{
					$rand_coordinates1[$suma] += $array1[$position_sum + $suma];#here we have to correct in case there are no chains in PDB
				}
				$number_atoms1 += 1;
			}
			if($aa_real == $aa_rand2)
			{
				for(my $suma = 0; $suma <= 2; $suma++)
				{
					$rand_coordinates2[$suma] += $array1[$position_sum + $suma];
				}
				$number_atoms2 += 1;
			}
			$linea += 1;
		}
		for(my $ql = 0; $ql <= 2; $ql++)
		{
			$rand_coordinates1[$ql] /= $number_atoms1;
			$rand_coordinates2[$ql] /= $number_atoms2;
		}			
		$rand_vector[$repeat] = Dimen::Atom_distance(\@rand_coordinates1, \@rand_coordinates2);
		$repeat += 1;
	}
	$Mean_atom_dist = Statistics::Average(@rand_vector);
	$SE_atom_dist = Statistics::SE(\@rand_vector, $Mean_atom_dist);
}##Here we still have to store the three-coordinates
my($n) = 0;

###STORING Z-TABLE VALUES INTO AN ARRAY
my($Z_file) = 'Z_scores';
my($z_man);
$z_man = File_reader::Input_file($Z_file);
my(@Z_scores);
my($score_value) = 0;
while(<$z_man>)
{
	my(@array) = split(/\t+/, $_);
	my($i) = 1;
	for(my($j) = 0; $j <= 9; $j++)
	{
		$Z_scores[$score_value] = $array[$i];
		$score_value += 1;
		$i += 1;
	}
}

###REAQDING SEQUENCES FROM INPUT FILE###
$information[0] =~ s/\t+//g;
print $out $header;
print $out "\nAmino acid sequence alignment for " . $information[0] . ":\n";
print $out "----------------------------------------------------------------------------\n";
#File_reader::Format_reader($information[4], $information[0], \@sequences, \@seq_names);
my($i) = 0;
while($i <= scalar(@sequences) - 1)
{
	print $out $seq_names[$i] . "\n". $sequences[$i] . "\n";
	$i += 1;
}
my($reference1, $reference2) = (0,0);
$reference1 = ($information[7] - 1);
if($information[3] == 1)
{
	$information[1] =~ s/\t+//g;
	#File_reader::Format_reader($information[5], $information[1], \@sequences2, \@seq_names2);
	$i = 0;
	print $out "\nAmino acid sequence alignment for " . $information[1] . ":\n";
	print $out "----------------------------------------------------------------------------\n";
	while($i <= scalar(@sequences) - 1)
	{
		print $out "$seq_names2[$i]\n$sequences2[$i] \n";
		$i += 1;
	}
	$reference2 = ($information[8] - 1);
}
if($reference1 >= scalar(@seq_names))
{
	print "Error: Could not find  $information[7] in the alignment\nPlease check case and control file!!!\n";
	exit;
}
if($reference2 >= scalar(@seq_names))
{
	print "Error: Could not find  $information[8] in the alignment\nPlease check case and control file!!!\n";
	exit;
}
File_reader::test_seqlength(@sequences);
###STORING BLOSUM VALUES###
my($BLOSUM) = 'Blosum';
my(@final_array_B, @Dif_per_seq, @Blosum_aa);
CAPS_module::Read_blosum(\$BLOSUM, \@final_array_B, \@Blosum_aa); #exhange between gaps and any other amino acid is marked by 0

###ESTIMATING PAIR-WISE AMINO ACID OR SYNONYMOUS NUCLEOTIDE DISTANCES###
my(@distances, @distances2, $distance, @distances_sorted, @distances_sorted2, @relative_dist, @relative_dist2, $relative);
if($information[13] == 1)
{
	my($loop1, $loop2) = (0,1);
	$i = 0;
	if($information[14] == 0)
	{
		print "\nEstimating Li-based pairwise synonymous distances....................\n";
		print $out "\n\nLi pair-wise synonymous distances\n";
		print $out "---------------------------------------------------------------\n\n";
	}elsif($information[14] == 1)
	{
		print "\nEstimating Poisson-corrected amino acid distances.......................\n";
		print $out "\n\nPair-wise Poisson-corrected amino acid distances\n";
		print $out "----------------------------------------------------------------------------------------\n\n";
	}	
	while($loop1 <= scalar(@seq_names) - 2)
	{
		print $out $seq_names[$loop1] . "  ";
		while($loop2 <= scalar(@seq_names) - 1)
		{
			if($information[14] == 0)
			{
				$distance = Distances::Li_synonymous($nuc_seq1[$loop1], $nuc_seq1[$loop2]);		
				if($information[3] == 1)
				{
					my($distance2) = Distances::Li_synonymous($nuc_seq2[$loop1], $nuc_seq2[$loop2]);
					$distances2[$i] = sprintf"%.4f", $distance2;
				}
				$distances[$i] = sprintf"%.4f", $distance;
			}elsif($information[14] == 1)
			{
				$distance = Distances::Poisson_dist($sequences[$loop1], $sequences[$loop2]);
				if($information[3] == 1)
				{
					my($distance2) = Distances::Poisson_dist($sequences2[$loop1], $sequences2[$loop2]);
					$distances2[$i] = sprintf"%.4f", $distance2;
				}
				$distances[$i] = sprintf"%.4f", $distance;
			}
			print $out $distances[$i] . " ";
			$i += 1;
			$loop2 += 1;
		}
		print $out "\n";
		$loop1 += 1;
		$loop2 = $loop1 + 1;
	}
	@distances_sorted = sort {$a <=> $b} @distances;
	my(@distances_sorted2);
	if($information[3] == 1)
	{
		@distances_sorted2 = sort {$a <=> $b} @distances2;
	}
	my($seqs) = 0;
	if($information[3] == 0)
	{
		Distances::Relative_distance(\@distances, \@distances_sorted, \@relative_dist);
	}else
	{
		Distances::Relative_distance(\@distances, \@distances_sorted, \@relative_dist);
		Distances::Relative_distance(\@distances2, \@distances_sorted2, \@relative_dist2);
	}
}
my(@aa_dif_column1, @aa_dif_column2, @aa_diferentes1, @aa_diferentes2);

###ESTIMATING INTRAMOLECULAR COEVOLUTION###
##In the co-evolution we need to store the information regarding Theta; Divergence corrected by the time;
my(@D, @D_corrected, @D_aver, @D_aver_red, @D_aver_red2, @D2, @D_corrected2, @D_aver2, @Correlations);
$n = 1;
my($m,  $number_corr) = (0,0);
print $out "\n\nCO-EVOLUTION ANALYSIS";
print $out "\n---------------------------\n\n";
if(($information[3] == 1) || (($information[3] == 0) && ($information[6] == 1)))
{
	print $out "Posicion AA1\tPosicion AA2\tMean D1\t\tMean D2\t\tCorrelation\n";
	print $out "------------\t------------\t-------\t\t-------\t\t-----------\n";
}
if(($information[3] == 0) &&  ($information[6] == 0))
{
	my($Mean_rand_red, $SE_rand_red);
	$Mean_rand_red = sprintf "%.4f", $Mean_atom_dist;
	$SE_rand_red = sprintf "%.4f", $SE_atom_dist;
	print $out "Mean random distance = $Mean_rand_red (SE = $SE_rand_red)\n\n";
	print $out "Posicion AA1\tPosicion AA2\tMean D1\t\tMean D2\t\tCorrelation\t\tAtom distance\n";
	print $out "------------\t------------\t-------\t\t-------\t\t-----------\t\t-------------\n";
}
my(@group, @number_aas, $D_mean, $D_mean2, $SE_D, $SE_D2);
if($information[3] == 0) #Case of intra-molecular co-evolution asked by the user
{
	my(@aa_differences, @groups_of_coevolution);
	Distances::Diff_aa_column(\@sequences, \@aa_differences);
	CAPS_module::Estimate_D(\@sequences,  \@D, \@D_corrected, \@relative_dist, \$information[13], \@final_array_B, \@Blosum_aa, \@aa_differences);
	Distances::Diff_aa_column(\@sequences, \@aa_dif_column1);##Calculating number of aa differences per column
	my($columna, $aa_diff) = (0,0);
	while($columna <= scalar(@aa_dif_column1) - 1)##storing columns with more than 1 different amino acids
	{
		if($aa_dif_column1[$columna] > 0)
		{
			$aa_diferentes1[$aa_diff] = $aa_dif_column1[$columna];
			$aa_diff += 1;
		}
		$columna += 1;
	}
	my($average_aa_diff) = Statistics::Average(@aa_diferentes1);
	$number_sequences = scalar(@sequences);
	my($length_seq) = length($sequences[0]);
	print "\nComputing intra-molecular co-evolution analysis, it might take few minutes, please wait...........\n";
	CAPS_module::Co_evolution_intra(\@D, \@D_corrected, \@Correlations, \$information[13], \$number_sequences, \$length_seq);
	if($information[13] == 0)
	{
		CAPS_module::D_average(\@D, \@D_aver, \$number_sequences);
		for(my $k = 0; $k <= scalar(@D_aver) - 1; $k++)
		{
			$D_aver_red[$k] = sprintf"%.4f", $D_aver[$k];
		}
	}elsif($information[13] == 1)
	{
		CAPS_module::D_average(\@D_corrected, \@D_aver, \$number_sequences);
		for(my $k = 0; $k <= scalar(@D_aver) - 1; $k++)
		{
			$D_aver_red[$k] = sprintf"%.4f", $D_aver[$k];
		}
	}
	my($cogroup) = 0;
	my(@coevolvings);
	my(@number_aas) = 0;
	my($threshold) = 0;
	##Deciding the method to detect significant co-evolution
	if($information[11] == 0)
	{
		$threshold = $information[12];
	}
	if($information[11] == 1)
	{
		CAPS_module::Get_threshold(\$information[15],\$information[16],\$threshold, \@Correlations);
	}else
	{
		#system '../Protseq-gen'; #Here we have to prepare all the information needed'
		#CAPS_module::Gets_means(simulated_file);
		#CAPS_module::Gets_threshold(\$threshold, \@simul_correlations);
	}
	my(@permut_corr, @permut_corr_ord);
	Get_permutational_F(\@Correlations, \@permut_corr);
	@permut_corr_ord = sort {$b <=> $a} @permut_corr;
	my(@ordered_correlations);
	@ordered_correlations = sort {$b <=> $a} @Correlations;
	$D_mean = Statistics::Average(@D_aver);# Is going to be used to select sites with sufficient signal for coevolution analysis
	$SE_D = Statistics::SE(\@D_aver, $D_mean);
	close $out;
	my($sizeGroup) = ($information[19]/100) * length($sequences[0]);
	print "\nGenerating the Groups of Coevolution.....................\n";
	if($information[6] == 0)
	{
		print_to_table(\@D, \@D2, \@D_aver_red, \@D_aver_red2, \@Correlations, \@coevolvings, \@number_aas, \$cogroup, \$threshold, \$information[2], \$D_mean, \$D_mean2, \@threeD_coordinates, \$sequences[$reference1], \$information[6], \$average_aa_diff, \@sequences, \$sizeGroup, \@groups_of_coevolution,\$information[18], \@ordered_correlations, \@permut_corr_ord);
	}else
	{
		my(@num) = 0;
		print_to_table(\@D, \@D2, \@D_aver_red, \@D_aver_red2, \@Correlations, \@coevolvings, \@number_aas, \$cogroup, \$threshold, \$information[2], \$D_mean, \$D_mean2, \@num, \$sequences[$reference1], \$information[6], \$average_aa_diff, \@sequences, \$sizeGroup, \@groups_of_coevolution,\$information[18], \@ordered_correlations, \@permut_corr_ord);
	}
	if(scalar(@groups_of_coevolution) > 0)
	{
		print "\nAnalysis of Compensatory Mutations:\nAnalysis of correlated Hydrophobicity.....\nAnalysis of correlated Molecular Weight.....\n";
		##DETECTION AND ANALYSIS OF COMPENSATORY MUTATIONS##
		my(@Hydros_per_site, @hydrophobicities, @MW);
		my($in_hydro) = "Hydrophobicity";
		open(H, $in_hydro) || die "Error: Could not find the file $in_hydro !!\n";
		@hydrophobicities = <H>;
		CompMut::Hydrophobicity(\@sequences, \@groups_of_coevolution, \@hydrophobicities, \@Hydros_per_site, \@MW);
		
		###Building a distribution of hydrophobic comparisons between pairs of sites###
		my(@groups_of_randoms, @Hydros_randoms, @MW_randoms);
		my($random) = 0;
		while($random <= 20000)
		{
			my($rand_number1, $rand_number2);
			$rand_number1 = int(rand(length($sequences[0])));
			do
			{
				$rand_number2 = int(rand(length($sequences[0])));
			}until($rand_number2 != $rand_number1);
			$groups_of_randoms[$random] = $rand_number1;
			$groups_of_randoms[$random + 1] = $rand_number2;
			$random += 2;
		}
		Calculate_hydros(\@sequences, \@groups_of_randoms, \@hydrophobicities, \@Hydros_randoms, \@MW_randoms);
		
		
		###TRANSFORMING HYDROPHOBICITY LEVELS INTO VARIABILITY OF HYDROPHOBICITY PER SITE IN EACH COEVOL. GROUP###
		my(@Hv, @Hv_randoms, @MWv, @MWv_randoms); #Hv = Hydrophobicity variability; MWv = Molecular weight variability
		my($num_group_coevol) = scalar(@groups_of_coevolution);
		my($num_seqs) = scalar(@sequences);
		my($num_groups_rand) = 10000;
		CompMut::HydroVar(\@Hydros_per_site, \$num_group_coevol, \$num_seqs, \@groups_of_coevolution, \@Hv, \@MW, \@MWv);
		Calculate_HydrosVar(\@groups_of_randoms, \@Hydros_randoms, \@Hv_randoms, \@sequences, \@MW_randoms, \@MWv_randoms);
		my(@HydCorr, @HydCorr_randoms, @MWCorr, @MWCorr_randoms);
		CompMut::HydroCorr(\@groups_of_coevolution, \@HydCorr, \@Hv, \@MWCorr, \@MWv);
		HydroCorr_randoms(\@groups_of_randoms, \@HydCorr_randoms, \@Hv_randoms, \@MWCorr_randoms, \@MWv_randoms);
		my(@sorted_Corr_rand) = sort {$a <=> $b} @HydCorr_randoms;
		my(@sorted_MWCorr_rand) = sort {$a <=> $b} @MWCorr_randoms;
		$i = 0;
		open(OUTPUT, ">>$information[2]");
		print OUTPUT "\n\tDetection of Compensatory Mutations\n\t------------------------------------------------------------\n\n";
		##PRINTING THE COEVOLVING PAIRS WITH CORRELATION IN THEIR HYDROPHOBICITIES##
		print OUTPUT "Detecting correlation in hydrophobicities in the groups of coevolution:\n\n";
		print OUTPUT "Group\t\t\tSites\t\tHydroph. Corr.\t\tProb\n";
		print OUTPUT "-----\t\t\t-----\t\t--------------\t\t----\n";
		my($Hydro_out) = "Hydros_coevol.csv";
		open(Hout, ">$Hydro_out");
		print Hout "Group,Site1,Site2,Hydro_Correl,Prob\n";
		my(@array_for_hydros);
		my($pointer_to_hydros) = 0;
		while($i <= scalar(@groups_of_coevolution) - 1)
		{
			my(@array_each_group) = split(/\s+/, $groups_of_coevolution[$i]);
			my($buc1, $buc2, $found) = (0,1,0);
			my($correlation_order) = 0;
			while($buc1 <= scalar(@array_each_group) - 2)
			{
				while($buc2 <= scalar(@array_each_group) - 1)
				{
					my($corr_prob, $corr_prob2) = (0,0);
					Locate_position(\@sorted_Corr_rand, \$HydCorr[$i][$correlation_order], \$corr_prob);
					my($corr_round, $prob_round);
					$corr_round = sprintf "%.4f", $HydCorr[$i][$correlation_order];
					$prob_round = sprintf "%.4f", $corr_prob;
					if($prob_round <= 0.05)
					{
						$found += 1;
						if($found >= 1)
						{
							print OUTPUT "G" . ($i + 1);
							print Hout ($i + 1) . ",";
						}
						print OUTPUT "\t\t" . $array_each_group[$buc1] . " & " . $array_each_group[$buc2] . "\t\t" . $corr_round . "\t\t$prob_round\n";
						print Hout $array_each_group[$buc1] . "," . $array_each_group[$buc2] . "," . $corr_round . "," . "$prob_round\n";
						$array_for_hydros[$pointer_to_hydros] = $array_each_group[$buc1] . " " . $array_each_group[$buc2];
						$pointer_to_hydros += 1;
					}
					$correlation_order += 1;
					$buc2 += 1;
				}
				$buc1 += 1;
				$buc2 = $buc1 + 1;
			}
			$i += 1;
			if($found > 0)
			{
				print OUTPUT "\n";
			}
		}
		##PRINTING THE COEVOLVING PAIRS WITH SIGNIFICANT CORRELATIONS IN THEIR MOLECULAR WEIGHTS##
		print OUTPUT "Detecting correlation in the Molecular weights in the groups of coevolution:\n\n";
		print OUTPUT "Group\t\t\tSites\t\tMW Corr.\tProb\n";
		print OUTPUT "-----\t\t\t-----\t\t--------\t-----\n";
		$i = 0;
		my($MW_out) = "MW_coevol.csv";
		open(MWout, ">$MW_out");
		print MWout "Group,Site1,Site2,MW_Correl,Prob\n";
		my(@array_for_MW);
		my($pointer_to_MW) = 0;
		while($i <= scalar(@groups_of_coevolution) - 1)
		{
			my(@array_each_group) = split(/\s+/, $groups_of_coevolution[$i]);
			my($buc1, $buc2, $found) = (0,1,0);
			my($correlation_order2) = 0;
			while($buc1 <= scalar(@array_each_group) - 2)
			{
				while($buc2 <= scalar(@array_each_group) - 1)
				{
					my($corr_prob2) = 0;
					Locate_position(\@sorted_MWCorr_rand, \$MWCorr[$i][$correlation_order2], \$corr_prob2);
					my($MWcorr_round, $MWprob_round);
					$MWcorr_round = sprintf "%.4f", $MWCorr[$i][$correlation_order2];
					$MWprob_round = sprintf "%.4f", $corr_prob2;
					if($MWprob_round <= 0.05)
					{
						$found += 1;
						if($found >= 1)
						{
							print OUTPUT "G" . ($i + 1);
							print MWout ($i + 1) . ",";
						}
						print OUTPUT "\t\t" . $array_each_group[$buc1] . " & " . $array_each_group[$buc2] . "\t" . $MWcorr_round . "\t\t$MWprob_round\n";
						print MWout $array_each_group[$buc1] . "," . $array_each_group[$buc2] . "," . $MWcorr_round . ",$MWprob_round\n";
						$array_for_MW[$pointer_to_MW] = $array_each_group[$buc1] . " " . $array_each_group[$buc2];
						$pointer_to_MW += 1;
					}
					$correlation_order2 += 1;
					$buc2 += 1;
				}
				$buc1 += 1;
				$buc2 = $buc1 + 1;
			}
			$i += 1;
			if($found > 0)
			{
				print OUTPUT "\n";
			}
		}
		
		##PRINTING THE COEVOLVING PAIRS WITH SIGNIFICANT CORRELATIONS IN THEIR MOLECULAR WEIGHTS AND HYDROPHOBICITIES##
		print OUTPUT "Detecting correlation in the Molecular weights and hydrophobicities in the groups of coevolution:\n\n";
		print OUTPUT "Group\t\t\tSites\t\tMW Corr.\tProb\t\tHydroph. Corr.\t\tProb\n";
		print OUTPUT "-----\t\t\t-----\t\t--------\t----\t--------------\t\t-----\n";
		$i = 0;
		my($HydMW_out) = "CompCoevol.csv";
		open(HMWout, ">$HydMW_out");
		print HMWout "\t\tSUMMARY OF CO-EVOLUTION ANALYSIS AND THE CORRELATED VARIATION OF CO-EVOLVING AMINO ACIDS IN THEIR HYDROPHOBICITY AND MOLECULAR WEIGHT\n\t\t1 indicates significant correlation; 0 indicates no significant correlation has been detected\n\n"; 
		print HMWout "Group,Site1,Site2,Coevolution,HydrCov,MWCov\n";
		while($i <= scalar(@groups_of_coevolution) - 1)
			{
			my(@array_each_group) = split(/\s+/, $groups_of_coevolution[$i]);
			my($buc1, $buc2, $found) = (0,1,0);
			my($correlation_order, $correlation_order2) = (0,0);
			while($buc1 <= scalar(@array_each_group) - 2)
			{
				while($buc2 <= scalar(@array_each_group) - 1)
				{
					my($corr_prob, $corr_prob2) = (0,0);
					Locate_position(\@sorted_Corr_rand, \$HydCorr[$i][$correlation_order], \$corr_prob);
					Locate_position(\@sorted_MWCorr_rand, \$MWCorr[$i][$correlation_order2], \$corr_prob2);
					my($corr_round, $prob_round, $MWcorr_round, $MWprob_round);
					$corr_round = sprintf "%.4f", $HydCorr[$i][$correlation_order];
					$prob_round = sprintf "%.4f", $corr_prob;
					$MWcorr_round = sprintf "%.4f", $MWCorr[$i][$correlation_order2];
					$MWprob_round = sprintf "%.4f", $corr_prob2;
					print HMWout ($i + 1) . "," . $array_each_group[$buc1] . "," . $array_each_group[$buc2] . ",1,";
					if(($prob_round <= 0.05) && ($MWprob_round <= 0.05))
					{
						$found += 1;
						if($found >= 1)
						{
							print OUTPUT "G" . ($i + 1);
						}
						print OUTPUT "\t\t" . $array_each_group[$buc1] . " & " . $array_each_group[$buc2] . "\t" . $MWcorr_round . "\t\t$MWprob_round\t\t $corr_round\t\t\t$prob_round\n";
					}
					my($search_in_array, $found_in_array) = (0,0);
					my($joined_sites) = $array_each_group[$buc1] . " " . $array_each_group[$buc2];
					while($search_in_array <= scalar(@array_for_hydros) - 1)
					{
						if($joined_sites eq $array_for_hydros[$search_in_array])
						{
							$found_in_array += 1;
						}
						$search_in_array += 1;
					}
					if($found_in_array > 0)
					{
						print HMWout "1,";
					}else
					{
						print HMWout "0,";
					}
					$search_in_array = 0; $found_in_array = 0;
					while($search_in_array <= scalar(@array_for_MW) - 1)
					{
						if($joined_sites eq $array_for_MW[$search_in_array])
						{
							$found_in_array += 1;
						}
						$search_in_array += 1;
					}
					if($found_in_array > 0)
					{
						print HMWout "1\n";
					}else
					{
						print HMWout "0\n";
					}
					$correlation_order += 1;
					$correlation_order2 += 1;
					$buc2 += 1;
				}
				$buc1 += 1;
				$buc2 = $buc1 + 1;
			}
			$i += 1;
			if($found > 0)
			{
				print OUTPUT "\n";
			}
		}
	}
close OUTPUT;
}else			##Case of Intermolecular co-evolution asked by the user
{
	my(@aa_differences1, @aa_differences2, @groups_of_coevolution);
	Distances::Diff_aa_column(\@sequences, \@aa_differences1);
	Distances::Diff_aa_column(\@sequences2, \@aa_differences2);
	CAPS_module::Estimate_D(\@sequences,  \@D, \@D_corrected, \@relative_dist, \$information[13], \@final_array_B, \@Blosum_aa, \@aa_differences1);
	CAPS_module::Estimate_D(\@sequences2,  \@D2, \@D_corrected2, \@relative_dist2, \$information[13], \@final_array_B, \@Blosum_aa,\@aa_differences2);
	###############
	Distances::Diff_aa_column(\@sequences, \@aa_dif_column1);##Calculating number of aa differences per column
	Distances::Diff_aa_column(\@sequences2, \@aa_dif_column2);
	my($columna, $aa_diff) = (0,0);
	while($columna <= scalar(@aa_dif_column1) - 1)##storing columns with more than 1 different amino acids
	{
		if($aa_dif_column1[$columna] > 0)
		{
			$aa_diferentes1[$aa_diff] = $aa_dif_column1[$columna];
			$aa_diff += 1;
		}
		$columna += 1;
	}
	my($average_aa_diff) = Statistics::Average(@aa_diferentes1);
	$columna = 0; $aa_diff = 0;
	while($columna <= scalar(@aa_dif_column2) - 1)##storing columns with more than 1 different amino acids
	{
		if($aa_dif_column2[$columna] > 0)
		{
			$aa_diferentes2[$aa_diff] = $aa_dif_column2[$columna];
			$aa_diff += 1;
		}
		$columna += 1;
	}
	my($average_aa_diff2) = Statistics::Average(@aa_diferentes2);
	$number_sequences = scalar(@sequences);
	print "\nComputing inter-molecular co-evolution analysis, it might take few minutes, please wait....................\n";
	CAPS_module::Co_evolution_inter(\@D, \@D2, \@D_corrected, \@D_corrected2, \@Correlations, \$information[13], \$number_sequences);
	my(@D_aver_red, @D_aver_red2);
	if($information[13] == 0)
	{
		CAPS_module::D_average(\@D, \@D_aver, \$number_sequences);
		for(my $k = 0; $k <= scalar(@D_aver) - 1; $k++)
		{
			$D_aver_red[$k] = sprintf"%.4f", $D_aver[$k];
		}
		CAPS_module::D_average(\@D2, \@D_aver2, \$number_sequences);
		for(my $k = 0; $k <= scalar(@D_aver2) - 1; $k++)
		{
			$D_aver_red2[$k] = sprintf"%.4f", $D_aver2[$k];
		}		
	}elsif($information[13] == 1)
	{
		CAPS_module::D_average(\@D_corrected, \@D_aver, \$number_sequences);
		for(my $k = 0; $k <= scalar(@D_aver) - 1; $k++)
		{
			$D_aver_red[$k] = sprintf"%.4f", $D_aver[$k];
		}
		CAPS_module::D_average(\@D_corrected2, \@D_aver2, \$number_sequences);
		for(my $k = 0; $k <= scalar(@D_aver2) - 1; $k++)
		{
			$D_aver_red2[$k] = sprintf"%.4f", $D_aver2[$k];
		}
	}
	my($cogroup) = 0;
	my(@coevolvings);
	my(@number_aas) = 0;
	my($threshold) = 0;
	##Deciding the method to detect significant co-evolution
	if($information[11] == 0)
	{
		$threshold = $information[12];
	}elsif($information[11] == 1)
	{
		CAPS_module::Get_threshold(\$information[15],\$information[16],\$threshold, \@Correlations);
	}else
	{
		#system '../Protseq-gen'; #Here we have to prepare all the information needed'
		#CAPS_module::Gets_means(simulated_file);
		#CAPS_module::Gets_threshold(\$threshold, \@simul_correlations);
	}
	$D_mean = Statistics::Average(@D_aver);# Is going to be used to select sites with sufficient signal for coevolution analysis
	$D_mean2 = Statistics::Average(@D_aver2);
	$SE_D = Statistics::Average(\@D_aver, $D_mean);
	$SE_D2 = Statistics::Average(\@D_aver2, $D_mean2);
	close $out;
	my($sizeGroup) = ($information[19]/100) * length($sequences[0]);
	print_table_inter(\@D, \@D2, \@D_aver_red, \@D_aver_red2, \@Correlations, \@coevolvings, \@number_aas, \$cogroup, \$threshold, \$information[2], \$D_mean, \$D_mean2, \$average_aa_diff, \$average_aa_diff2, \@sequences, \@sequences2, \$sequences[$reference1],  \$sequences2[$reference2], \$information[18], \$sizeGroup, \@groups_of_coevolution);
	
}		
exit;

###SUBROUTINES USED IN THE MAIN BODY OF THE PROGRAM

#SUBROUTINE TO PRINT THE DATA INTO THE TABLE
sub print_to_table
{
	my($D, $D2, $D_aver_red, $D_aver_red2, $Correlations, $coevolvings, $number_aas, $cogroup, $information, $information_out, $MEAN, $MEAN2, $threeD_coor, $sequence_reference, $information_test, $average_difference, $sequences, $sizeG, $G_of_coevolution, $Min_R, $ord_corr, $permutations) = @_;
	my($out);
	open($out, ">>$$information_out");
	my($n, $number_corr, $atom_distance) = (1,0,0);
	my(@coevolvings_real, @pairs_of_sites);
	my($pairs_pointer) = 0;
	my($percent_comp, $point_to_percent) = (0,0);
	while($m <= scalar(@$D) - 2)
	{
		my($aa_diff_pos1, $aa_diff_pos2, $aa1_pars, $aa2_pars, $aa1_parsimonity) = (0,0,0,0,0);
		my($seq_number) = 0;
		my(@sequences_pos, @diferencias);
		##Estimating the number of differences per column
		while($seq_number <= scalar(@$sequences) - 1)
		{
			 $sequences_pos[$seq_number] = substr($$sequences[$seq_number], $m, 1);
			 $seq_number += 1;
		}
		Distances::Diff_aa_column(\@sequences_pos, \@diferencias);
		my(@seqs);
		my($point_to_seqs) = 0;
		while($point_to_seqs < scalar(@$sequences))
		{
			$seqs[$point_to_seqs] = $$sequences[$point_to_seqs];
			$point_to_seqs += 1;
		}
		$aa_diff_pos1 = $diferencias[0];
		##################################
		##CALCULATING THE REAL POSITION IN THE REFERENCE SEQUENCE##
		my($pos1) = $m + 1;
		Parsimonity(\$m, \@seqs, \$aa1_parsimonity);
		my($l) = 0;
		$$coevolvings[$$cogroup][$l] = $pos1;
		
		my($pos, $gaps1) = (0,0);
		while($pos <= $pos1 - 1)
		{
			my($aa) = substr($$sequence_reference, $pos, 1);
			if($aa eq '-')
			{
				$gaps1 += 1;
			}
			$pos += 1;
		}
		my($final_pos1) = $pos1 - $gaps1;
		$coevolvings_real[$$cogroup][$l] = $final_pos1;
		$l += 1;
		$$number_aas[$$cogroup] = 0;
		$$number_aas[$$cogroup] += 1;
		my($final_coevol_pointer) = 0;
		my(@final_coevolvings, @Super_final_coevolvings);
		while($n <= scalar(@$D) - 1)
		{
			$percent_comp  = (($point_to_percent)/((scalar(@$D)) * (scalar(@$D) - 1)/2)) * 100;
			my($per_red) = sprintf "%.0f", $percent_comp;
			print "\r$per_red%";
			$seq_number = 0;
			my($aa2_parsimonity) = 0;
			##Estimating the number of differences per column
			while($seq_number <= scalar(@$sequences) - 1)
			{
				 $sequences_pos[$seq_number] = substr($$sequences[$seq_number], $n, 1);
				 $seq_number += 1;
			}
			Distances::Diff_aa_column(\@sequences_pos, \@diferencias);
			$aa_diff_pos2 = $diferencias[0];
			####################################
			my($Correlation_red);
			$Correlation_red = sprintf"%.4f", $$Correlations[$number_corr];
			my($pos2) = $n + 1;
			Parsimonity(\$n, \@seqs, \$aa2_parsimonity);
			#At this point we should compute the significance of D values with a subroutine
			my($gaps1, $gaps2) = (0,0);
			my(@coordinates1, @coordinates2);
			$pos = 0;
			while($pos <= $pos2 - 1)
			{
				my($aa) = substr($$sequence_reference, $pos, 1);
				if($aa eq '-')
				{
					$gaps2 += 1;
				}
				$pos += 1;
			}
			my($final_pos2) = $pos2 - $gaps2;
			if($final_pos1 == $final_pos2)
			{
				$final_pos2 = $final_pos1 + 1; ##This is the case where we have gaps in the reference sequence
			}
			#if((($Correlation_red >= $$information) || ($Correlation_red * (-1) >= $$information)) && ($$D_aver_red[$m] > $$MEAN) && ($$D_aver_red[$n] > $$MEAN) && ($aa_diff_pos1 > $$average_difference) && ($aa_diff_pos2 > $$average_difference))
			#my($limit) = int(scalar(@$sequences) * 0.1);
			my($pos_corr, $ql) = (0,0);
			if($Correlation_red >= $$information)
			{
				while($ql <= scalar(@$permutations) - 1)
				{
					if($$Correlations[$number_corr] >= $$ord_corr[$pos_corr])
					{
						last;
					}else
					{
						$pos_corr += 1;
					}
					$ql += 1;
				}
			}
			if($pos_corr <= scalar(@$permutations) - 1)
			{
				if($Correlation_red < $$permutations[$pos_corr])
				{
					$Correlation_red = 0;
				}
			}
			if((($Correlation_red >= $$information) || ($Correlation_red * (-1) >= $$information)) && ($aa_diff_pos1 >= 2) && ($aa_diff_pos2 >= 2) && ($Correlation_red >= $$Min_R) && ($aa1_parsimonity == 1) && ($aa2_parsimonity == 1) )
			{
				$final_coevolvings[$final_coevol_pointer] = '';
				if($$information_test == 0)
				{
					my(@array_coor, @array_coor1, @array_coor2);
					my($search, $found) = (0,0);
					while($search <= scalar(@$threeD_coor) - 1)
					{
						@array_coor = split(/\s+/, $$threeD_coor[$search]);
						if($array_coor[0] == $final_pos1)
						{
							$array_coor1[0] = $array_coor[1]; $array_coor1[1] = $array_coor[2]; $array_coor1[2] = $array_coor[3];
							$found += 1;
						}
						if($array_coor[0] == $final_pos2)
						{
							$array_coor2[0] = $array_coor[1]; $array_coor2[1] = $array_coor[2]; $array_coor2[2] = $array_coor[3];
							$found += 1;
					#		last;
						}
						if($found == 2)
						{
							last;
						}
						$search += 1;
					}
					if($found > 1)
					{
						$atom_distance = Dimen::Atom_distance(\@array_coor1, \@array_coor2);
					}else
					{
						$atom_distance = 9999;
					}
					my($atom_dist_red) = sprintf "%.4f", $atom_distance;
					#print $out $pos1 . "(" . $final_pos1 . ")" . "\t\t" . $pos2 . "(" . $final_pos2 . ")" . "\t\t" . $$D_aver_red[$m] . "\t\t" . $$D_aver_red[$n] . "\t\t" . $Correlation_red . "\t\t\t" . $atom_dist_red . "\n";
					
					#############
					$final_coevolvings[$final_coevol_pointer] .= $pos1 . "(" . $final_pos1 . ")" . "\t" . $pos2 . "(" . $final_pos2 . ")" . "\t" . $$D_aver_red[$m] . "\t\t" . $$D_aver_red[$n] . "\t\t" . $Correlation_red . "\t\t\t" . $atom_dist_red;
					$final_coevol_pointer += 1;
					##############
					
							
				}else
				{
			#		print $out $pos1 . "(" . $final_pos1 . ")" . "\t\t" . $pos2 . "(" . $final_pos2 . ")" . "\t\t" . $$D_aver_red[$m] . "\t\t" . $$D_aver_red[$n] . "\t\t" . $Correlation_red . "\n";
					
					#############
					$final_coevolvings[$final_coevol_pointer] .= $pos1 . "(" . $final_pos1 . ")" . "\t" . $pos2 . "(" . $final_pos2 . ")" . "\t" . $$D_aver_red[$m] . "\t\t" . $$D_aver_red[$n] . "\t\t" . $Correlation_red;
					$final_coevol_pointer += 1;
					###############
					
				}
				$$coevolvings[$$cogroup][$l] = $pos2;
				$coevolvings_real[$$cogroup][$l] = $final_pos2;
				$$number_aas[$$cogroup] += 1;
				$l += 1;
			}
			$number_corr += 1;
			$point_to_percent += 1;
			$n += 1;
		}
		if((scalar(@final_coevolvings) <= $$sizeG) && (scalar(@final_coevolvings) > 0)) 
		{
			for(my $qul = 0; $qul <= scalar(@final_coevolvings) - 2; $qul++)
			{
				print $out $final_coevolvings[$qul] . "\n";
				my(@line) = split(/\t+/, $final_coevolvings[$qul]);
				$pairs_of_sites[$pairs_pointer] = "$line[0] $line[1]";
				$pairs_pointer += 1;
			}
		}
		$m += 1;
		$n = $m + 1;
		$$cogroup += 1;
	}
	if(scalar(@pairs_of_sites) > 0)
	{
		print $out "\nGroups of Co-evolving amino acids:";
		print $out "\n-------------------------------------------------------------\n";
	}
	my($groups, $group_pointer, $final_pointer, $pointer_to_final_group) = (0,0,0,0);
	my(@final_group, @subgroup1, @subgroup2);
	$final_group[0] = $pairs_of_sites[0];
	$groups += 1;
	while($groups <= scalar(@pairs_of_sites) - 1)
	{
		$final_pointer = 0;
		@subgroup2 = split(/\s/, $pairs_of_sites[$groups]);
		my($found1, $found2, $superfound) = (0,0,0);
		while($final_pointer <= $pointer_to_final_group)
		{
			@subgroup1 = split(/\s+/, $final_group[$pointer_to_final_group]);
			my($site1_pointer) = 0;
			$found1 = 0; $found2 = 0;
			$superfound = 0;
			while($site1_pointer <= scalar(@subgroup1) - 1)
			{
				if($subgroup2[0] eq $subgroup1[$site1_pointer])
				{
					$found1 += 1;
				}
				if($subgroup2[1] eq $subgroup1[$site1_pointer])
				{
					$found2 += 1;
				}
				$site1_pointer += 1;
			}
			if(($found1 > 0) && ($found2 == 0))
			{
				
				my($subdivisor, $found_at_last)  = (0,0);
				while($subdivisor <= scalar(@subgroup1) - 1)
				{
					my($pair) = $subgroup1[$subdivisor] . " " . $subgroup2[1];
					#my($busca) = $groups + 1;
					my($busca) = 0;#####
					
					while($busca <= scalar(@pairs_of_sites) - 1)
					{
						if($pairs_of_sites[$busca] eq $pair)
						{
							$found_at_last += 1;
							last;
						}
						$busca += 1;
					}
					#if($found_at_last > 0)
					#{
					#	$final_group[$pointer_to_final_group] .= " $subgroup2[1]";
					#	$superfound += 1;
					#	last;
					#}
					$subdivisor += 1;
				}
				#if($found_at_last > 0)
				#{
			#		last;
			#	}
				if($found_at_last >= scalar(@subgroup1))
				{
					$final_group[$pointer_to_final_group] .= " $subgroup2[1]";
					$superfound += 1;
				}
			
			}
			if(($found1 > 0) && ($found2 > 0))
			{
				$superfound += 1;
				last;
			}
			$final_pointer += 1;
		}	
		if($superfound == 0)
		{
		#	print $final_group[$pointer_to_final_group] . "\n";
			$pointer_to_final_group += 1;
			$final_group[$pointer_to_final_group] = $pairs_of_sites[$groups];
		}
		$groups += 1;	
	}
##Verifying and collapsing groups of coevolution####
	my(@G_final, @G1, @G2);
	$i = 1;
	my($pointer_to_G) = 1;
	$G_final[0] = $final_group[0];
	my($never_found) = 0;
	while($i <= scalar(@final_group) - 1)
	{
		@G1 = split(/\s+/, $final_group[$i]);
		my($j) = 0;
		my($new_group, $found_it) = (0,1);
		while($j <= scalar(@G_final) - 1)
		{
			@G2 = split(/\s+/, $G_final[$j]);
			my($k, $found_in) = (0,0);
			my(@differences);
			Find_G1_in_G2(\@G1, \@G2, \$found_in);#Mira si todos los elementos de G1 coinciden con los de G2
			if($found_in == scalar(@G1))##Si asi fuera
			{
				$found_it = 0;
				last;##Acaba el bucle
			}else
			{
				Build_array_differences(\@G1, \@G2, \@differences); ##Construye una array con las diferencias entre G1 y G2
				my($found_pair) = 0;
				
				
				Find_the_pairs(\@G2, \@differences, \@pairs_of_sites, \$found_pair);##Construye pares con cada elemento de G1 y de G2
				if($found_pair == scalar(@differences)) #Si encuentras los pares
				{
					Insert_difference_in_array(\$G_final[$j], \@differences);##anyades los elementos al grupo existente
					$found_it = 0;
					last;
				}
				$j += 1;
			}
		}
		if($found_it > 0)
		{
			$G_final[$pointer_to_G] = $final_group[$i];
			$pointer_to_G += 1;
		}
		$never_found = 0;
		$i += 1;
	}	
	$groups = 0;
	if((scalar(@G_final) > 0) && ($G_final[$groups]) && (length($G_final[$groups]) > 0))
	{
		while($groups <= scalar(@G_final) - 1)
		{
			print $out "G" . ($groups + 1) . ": " . $G_final[$groups] . "\n";
			$$G_of_coevolution[$groups] = $G_final[$groups];
			$groups += 1;
		}
	}else
	{
		print "There are no groups of coevolution\n\n";
	}		
}
###Table of significant correlations for inter-coevolution analysis
sub print_table_inter
{
	my($D, $D2, $D_aver_red, $D_aver_red2, $Correlations, $coevolvings, $number_aas, $cogroup, $information, $information_out, $MEAN, $MEAN2, $average_difference, $average_difference2, $sequences, $sequences2, $sequence_reference1, $sequence_reference2, $information_R, $sizeG, $G_of_coevolution) = @_;
	my($out);
	my($Coevolve_group) = 0;
	##Testing the values for the parameters
	my($test_out) = "Test.out";
	open(TEST, ">$test_out");
	print TEST "MEAN = $$MEAN; Mean2 = $$MEAN2;\tAverage_diff = $$average_difference; Average_diff2 = $$average_difference2\n\n";
	my(@coevolvings_real, @pairs_of_sites);
	my($pairs_pointer) = 0;
	
	open($out, ">>$$information_out");
	my($n, $number_corr, $m) = (0,0,0);
	while($m <= scalar(@$D) - 1)
	{
		my($aa_diff_pos1, $aa_diff_pos2) = (0,0);
		my(@sequences_pos, @diferencias);
		my($seq_number) = 0;
		##Estimating the number of differences per column
		while($seq_number <= scalar(@$sequences) - 1)
		{
			 $sequences_pos[$seq_number] = substr($$sequences[$seq_number], $m, 1);
			 $seq_number += 1;
		}
		Distances::Diff_aa_column(\@sequences_pos, \@diferencias);
		$aa_diff_pos1 = $diferencias[0];
		##################################
	
		my($pos1) = $m + 1;
		my($l) = 0;
		my($pos, $gaps1) = (0,0);
		while($pos <= $pos1 - 1)
		{
			my($aa) = substr($$sequence_reference1, $pos, 1);
			if($aa eq '-')
			{
				$gaps1 += 1;
			}
			$pos += 1;
		}
		my($final_pos1) = $pos1 - $gaps1;
		#$$coevolvings[$$cogroup][$l] = $pos1;
		$l += 1;
		#$$number_aas[$$cogroup] = 0;
		#$$number_aas[$$cogroup] += 1;
		my($final_coevol_pointer) = 0;
		my(@final_coevolvings);
		while($n <= scalar(@$D2) - 1)
		{
			$seq_number = 0;
			##Estimating the number of differences per column
			while($seq_number <= scalar(@$sequences2) - 1)
			{
				 $sequences_pos[$seq_number] = substr($$sequences2[$seq_number], $n, 1);
				 $seq_number += 1;
			}
			Distances::Diff_aa_column(\@sequences_pos, \@diferencias);
			$aa_diff_pos2 = $diferencias[0];
			####################################
			
			my($Correlation_red);
			$Correlation_red = sprintf"%.4f", $$Correlations[$number_corr];
			my($pos2) = $n + 1;
			$pos = 0;
			my($gaps2) = 0;
			while($pos <= $pos2 - 1)
			{
				my($aa) = substr($$sequence_reference2, $pos, 1);
				if($aa eq '-')
				{
					$gaps2 += 1;
				}
				$pos += 1;
			}
			my($final_pos2) = $pos2 - $gaps2;
			#if($final_pos1 == $final_pos2)
		#	{
		#		$final_pos2 = $final_pos1 + 1; ##This is the case where we have gaps in the reference sequence
		#	}
		
			#At this point we should compute the significance of D values with a subroutine
			#if((($Correlation_red >= $$information) ||  ($Correlation_red * (-1) >= $$information)) && ($$D_aver_red[$m] > $$MEAN) && ($$D_aver_red2[$n] > $$MEAN2) && ($aa_diff_pos1 >= $$average_difference) && ($aa_diff_pos2 >= $$average_difference2))
			#if((($Correlation_red >= $$information) ||  ($Correlation_red * (-1) >= $$information)) && ($aa_diff_pos1 >= $$average_difference) && ($aa_diff_pos2 >= $$average_difference2))
			#if((($Correlation_red >= $$information) ||  ($Correlation_red * (-1) >= $$information)) && ($aa_diff_pos1 >= 2) && ($aa_diff_pos2 >= 2) && (($Correlation_red >= $$information_R) ||  ($Correlation_red * (-1) >= $$information_R)))
			#if(($Correlation_red > 0) ||  ($Correlation_red * (-1) > 0)) 
			#my($limit) = int(scalar(@$sequences) * 0.1);
			if((($Correlation_red >= $$information) || ($Correlation_red * (-1) >= $$information)) && ($aa_diff_pos1 >= 2) && ($aa_diff_pos2 >= 2))
			{
				$final_coevolvings[$final_coevol_pointer] = '';
				print $out $pos1 . "($final_pos1)\t\t" . $pos2 . "($final_pos2)\t\t" . $$D_aver_red[$m] . "\t\t" . $$D_aver_red2[$n] . "\t\t" . $Correlation_red . "\n";
				$$coevolvings[$Coevolve_group][$l] = $pos2;
				$$number_aas[$Coevolve_group] += 1;
				$final_coevolvings[$final_coevol_pointer] .= $pos1 . "(" . $final_pos1 . ")" . "\t\t" . $pos2 . "(" . $final_pos2 . ")" . "\t\t" . $$D_aver_red[$m] . "\t\t" . $$D_aver_red[$n] . "\t\t" . $Correlation_red;
				$final_coevol_pointer += 1;
				$l += 1;
			}
			$number_corr += 1;
			$n += 1;
		}
	#	if($$number_aas[$Coevolve_group] > 0)
	#	{
	#		$Coevolve_group += 1;
	#	}
		if((scalar(@final_coevolvings) <= $$sizeG) && (scalar(@final_coevolvings) > 0)) 
		{
			for(my $qul = 0; $qul <= scalar(@final_coevolvings) - 2; $qul++)
			{
				print $out $final_coevolvings[$qul] . "\n";
				my(@line) = split(/\t+/, $final_coevolvings[$qul]);
				$pairs_of_sites[$pairs_pointer] = "$line[0] $line[1]";
				$pairs_pointer += 1;
			}
		}
		$m += 1;
		$n = 0;
		#$$cogroup += 1;
	}
	
	##BUILDING THE GROUPS OF COEVOLUTION##
	print $out "\nGroups of inter-molecular Co-evolving amino acids:";
	print $out "\n------------------------------------------------------------------------------------------\n";
	my($groups, $group_pointer, $final_pointer, $pointer_to_final_group, $number_of_pairs, $pairs, $subgroup_pointer) = (1,0,0,0,0,0,0);
	my(@sub_final_groupprot1, @sub_final_groupprot2, @final_group, @subgroup1, @subgroup2);
	$final_group[0] = $pairs_of_sites[0];
	my(@first_array) = split(/\s/, $pairs_of_sites[0]);
	$sub_final_groupprot1[0] = $first_array[0]; #This is going to build the groups in each protein
	$sub_final_groupprot2[0] = $first_array[1];
	$number_of_pairs = scalar(@pairs_of_sites);
	my(@array_for_groupings);
	while($pairs < $number_of_pairs)
	{
		@first_array = split(/\s/, $pairs_of_sites[$pairs]);
		$sub_final_groupprot1[$subgroup_pointer] = $first_array[0];
		$sub_final_groupprot2[$subgroup_pointer] = $first_array[1];
		my(@arrayA) = split(/\s+/, $sub_final_groupprot1[$subgroup_pointer]);
		$groups = 1 + $pairs;
		while($groups <= scalar(@pairs_of_sites) - 1)
		{
			my(@arrayB) = split(/\s/, $pairs_of_sites[$groups]);
			if($arrayA[0] eq $arrayB[0])
			{
				$sub_final_groupprot2[$subgroup_pointer] .= " " . $arrayB[1];
				$pairs += 1;
			}
			$groups += 1;
		}

		$array_for_groupings[$subgroup_pointer] = $sub_final_groupprot1[$subgroup_pointer] . " " . $sub_final_groupprot2[$subgroup_pointer];
		$subgroup_pointer += 1;
		$pairs += 1;
	}
 	my(%Coevol_groups);
	my($Gsize) = $$sizeG;
	%Coevol_groups = InterSort::make_final_groups(\@array_for_groupings, \$Gsize);
 	foreach my $element(keys %Coevol_groups)
	{
		print $out $element . "[" . $Coevol_groups{$element} . "]\n";
	}
 
 ##############################################################
 ##############################################################
  ##Here we are going to create the final groups of coevolution
  $subgroup_pointer = 0; $pairs = 0;
  my($subgroup_pointer2, $group_formed, $pointer_final_group) = (1,0,0);
  my(@final_group1, @final_group2, @temporal_group1, @temporal_group2);
  
}
##SUBROUTINE TO DEFINE THE MAIN GROUPS OF CORRELATED SITES
sub Select_groups
{
	my($selected_groups, $final_group) = @_;
	my($point_to_group) = 1;
	$$selected_groups[0] = $$final_group[0];
	my($i, $j) = (1,0);
	while($i <= scalar(@$final_group) - 1)
	{
		my(@array1, @array2, @location);
		@array1 = split(/\s+/, $$final_group[$i]);
		$j = 0;
		my($search, $found) = (0,0);
		while($j <= scalar(@$selected_groups) - 1)
		{
			@array2 = split(/\s+/, $$selected_groups[$j]);
			$location[$j] = 0;
			$search = 0;
			while($search <= scalar(@array2) - 1)
			{
				if($array1[0] == $array2[$search])
				{
					$found += 1;
				}
				if($found > 0)
				{
					$location[$j] += 1;
				}
				$search += 1;
			}
			$j += 1;
		}
		my($m, $n) = (0,0);
		while($m <= scalar(@location) - 1)
		{
			my(@array_fin) = split(/\s+/, $$final_group[$i]);
			if(scalar(@array_fin) <= $location[$m])
			{
				$n += 1;
			}
			$m += 1;
		}
		if(($found == 0) || ($n == 0))
		{
			$$selected_groups[$point_to_group] = $$final_group[$i];
			$point_to_group += 1;
		}
		$i += 1;
	}
}

##Subroutine to calculate hydrophobicities for random pairs
sub Calculate_hydros
{
	my($sequences, $groups_of_randoms, $hydrophobicities, $Hydros_randoms, $MW_randoms) = @_;
	my($i) = 0;
	while($i <= 19999)
	{
		my($seq_order) = 0;
		$$Hydros_randoms[$i] = '';
		$$MW_randoms[$i] = '';
		while($seq_order <= scalar(@$sequences) - 1)
		{
			my($aa);
			$aa = substr($$sequences[$seq_order], $$groups_of_randoms[$i], 1);
			my($search) = 0;
			while($search <= scalar(@$hydrophobicities) - 1)
			{
				my(@array2);
				@array2 = split(/\t+/, $$hydrophobicities[$search]);
				if($aa eq $array2[1])
				{
					if(length($$Hydros_randoms[$i]) > 0)
					{
						$$Hydros_randoms[$i] .= ' '  . $array2[2];
						$$MW_randoms[$i] .= ' ' . $array2[4];
					}else
					{
						$$Hydros_randoms[$i] = $array2[2];
						$$MW_randoms[$i] = $array2[4];
					}
					last;
				}
				if($aa eq "-")
				{
					if(length($$Hydros_randoms[$i]) > 0)
					{
						$$Hydros_randoms[$i] .= ' '  . 0;
						$$MW_randoms[$i] .= ' ' . 0;
					}else
					{
						$$Hydros_randoms[$i] = 0;
						$$MW_randoms[$i] = 0;
					}
					last;
				}
				$search += 1;
			}
			$seq_order += 1;
		}
		$i += 1;
	}
}
##Subroutine to calculate the variability in hydrophobicity in each site
sub Calculate_HydrosVar
{
	my($groups_of_randoms, $Hydros_randoms,  $Hv_randoms, $sequences, $MW_randoms, $MWv_randoms) = @_;
	my($i) = 0;
	while($i <= 19999)
	{
		my($buc1, $buc2) = (0,1);
		$$Hv_randoms[$i] = '';
		$$MWv_randoms[$i] = '';
		my(@array) = split(/\s+/, $$Hydros_randoms[$i]);
		my(@array_MW) = split(/\s+/, $$MW_randoms[$i]);
		while($buc1 <= scalar(@array) - 2)##was scalar(@$sequences) - 2
		{
			while($buc2 <= scalar(@array) - 1)##was scalar(@$sequences) - 1
			{
				if(length($$Hv_randoms[$i]) == 0)
				{
					$$Hv_randoms[$i] = sqrt(($array[$buc1] - $array[$buc2])**2);
					$$MWv_randoms[$i] = sqrt(($array_MW[$buc1] - $array_MW[$buc2])**2);
				}else
				{
					$$Hv_randoms[$i] .= " " . sqrt(($array[$buc1] - $array[$buc2])**2);
					$$MWv_randoms[$i] .= " " . sqrt(($array_MW[$buc1] - $array_MW[$buc2])**2);
				}
				$buc2 += 1;
			}
			$buc1 += 1;
			$buc2 = $buc1 + 1;
		}
		$i += 1;
	}
}
##Subroutine to calculate the correlation coeficient of 
sub HydroCorr_randoms
{
	my($groups_of_randoms, $HydCorr_randoms, $Hv_randoms, $MWCorr_randoms, $MWv_randoms) = @_;
	my($i, $j) = (0,0);
	while($i <= 19999)
	{
		my(@array1) = split(/\s+/, $$Hv_randoms[$i]);
		my(@array2) = split(/\s+/, $$Hv_randoms[$i + 1]);
		my(@arrayMW1) = split(/\s+/, $$MWv_randoms[$i]);
		my(@arrayMW2) = split(/\s+/, $$MWv_randoms[$i + 1]);
		$$HydCorr_randoms[$j] = Statistics::Correlation(\@array1, \@array2);
		$$MWCorr_randoms[$j] = Statistics::Correlation(\@arrayMW1, \@arrayMW2);
		if($$HydCorr_randoms[$j] < 0)
		{
			$$HydCorr_randoms[$j] *= (-1);
		}
		if($$MWCorr_randoms[$j] < 0)
		{
			$$MWCorr_randoms[$j] *= (-1);
		}
		$j += 1;
		$i += 2;
	}
}
##Subroutine to place a correlation in real position and calculate probability
sub Locate_position
{
	my($sorted_Corr_rand, $HydCorr, $prob) = @_;
	my($search) = 0;
	my($corr);
	if($$HydCorr < 0)
	{
		$corr = $$HydCorr * (-1);
	}else
	{
		$corr = $$HydCorr;
	}
	while($search <= scalar(@$sorted_Corr_rand) - 1)
	{
		if($corr <= $$sorted_Corr_rand[$search])
		{
			last;
		}
		$search += 1;
	}
	$$prob = (1 - $search/(scalar@$sorted_Corr_rand));
}
sub Find_G1_in_G2
{
	my($G1, $G2, $found_in) = @_;
	my($q) = 0;
	while($q <= scalar(@$G1) - 1)
	{
		my($K) = 0;
		while($K <= scalar(@$G2) - 1)
		{
			if($$G1[$q] eq $$G2[$K])
			{
				$$found_in += 1;
				last;
			}
			$K += 1;
		}
		$q += 1;
	}
}
sub Build_array_differences
{
	my($G1, $G2, $differences) = @_;
	my($q, $new_diff, $found, $dif) = (0,0,0,0);
	while($q <= scalar(@$G1) - 1)
	{
		my($K) = 0;
		while($K <= scalar(@$G2) - 1)
		{
			if($$G1[$q] eq $$G2[$K])
			{
				$found += 1;
				last;
			}
			$K += 1;
		}
		if($found == 0)
		{
			$$differences[$dif] = $$G1[$q];
			$dif += 1;
		}
		$q += 1;
		$found = 0;
	}
}
sub Find_the_pairs
{
	my($G2, $differences, $pairs_of_sites, $found_pair) = @_;
	my($K, $q) = (0,0);
	while($K <= scalar(@$differences) - 1)
	{
		$q = 0;
		my($found, $equivalente) = (0,0);
		while($q <= scalar(@$G2) - 1)
		{
			my($actual_pair);
			my($pairs) = 0;
			while($pairs <= scalar(@$pairs_of_sites) - 1)
			{
				my($search1, $search2) = (0,0);
				my(@array);
				@array = split(/\s+/, $$pairs_of_sites[$pairs]);
				if($$G2[$q] ne $$differences[$K])
				{
					if(($$G2[$q] eq $array[0]) || ($$G2[$q] eq $array[1]))
					{
						$search1 += 1;
					}
					if(($$differences[$K] eq $array[0]) || ($$differences[$K] eq $array[1]))
					{
						$search2 += 1;
					}
					if(($search1 > 0 ) && ($search2 > 0))
					{
						$found += 1;
						last;
					}
				}else
				{
					$equivalente += 1;
				}
				
				$pairs += 1;
			}
			$q += 1;
		}
		if($found >= (scalar(@$G2) - $equivalente))
		{
			$$found_pair += 1;
		}
		$K += 1;
	}
}
sub Insert_difference_in_array
{
	my($G_final, $differences) = @_;
	my($q, $found) = (0,0);
	while($q <= scalar(@$differences) - 1)
	{
		my($g) = 0;
		my(@array) = split(/\s+/, $$G_final);
		$found = 0;
		while($g <= scalar(@array) - 1)
		{
			if($$differences[$q] eq $array[$g])
			{
				$found += 1;
				last;
			}
			$g += 1;
		}
		if($found == 0)
		{
			my($G_group) = $$G_final;
			my($dif_to_insert) = $$differences[$q];
			Insert_in_order(\$G_group, \$dif_to_insert);
			$$G_final = $G_group;
		#	$$G_final .= " " . $$differences[$q];
		}
		$q += 1;
	}
}
sub Insert_in_order
{
	my($g, $dif) = @_;
	my(@array_g) = split(/\s+/, $$g);
	my($element) = 0;
	my(@array_dif) = split(/\(/, $$dif);
	my($inserted) = 0;
	my(@array);
	while($element <= scalar(@array_g) - 1)
	{
		@array = split(/\(/, $array_g[$element]);
		if(($array_dif[0] < $array[0]) && ($inserted == 0))
		{
			if($element == 0)
			{
				$$g = $$dif;
			}
			if($element > 0)
			{
				$$g .= " " . $$dif;
			}
			$inserted += 1;
		}
		if(($inserted > 0) || ($array_dif[0] > $array[0]))
		{
			if($element == 0)
			{
				$$g = $array_g[$element];
			}else
			{
				$$g .= " " . $array_g[$element];
			}
		}
		$element += 1;
	}
	if(($array_dif[0] > $array[0]) && ($inserted == 0))
	{
		$$g .= " " . $$dif;
	}
}
sub Parsimonity
{
	my($position, $seqs, $parsimonity) = @_;
	my($i) = 0;
	my(@aa);
	while($i <= scalar(@$seqs) - 1)
	{
		$aa[$i] = substr($$seqs[$i], $m, 1);
		$i += 1;
	}
	my($aa_diff) = 0;
	evaluate_differences(\@aa, \$aa_diff);
	if($aa_diff > 0)
	{
		$$parsimonity = 1;
	}
}	

sub evaluate_differences
{
	my($aa, $dif) = @_;
	my($k, $j) = (0,0);
	my(@aa_differences);
	$aa_differences[0] = $$aa[0];
	my($pointer_to_dif) = 1;
	my(@repetition);
	while($k <= scalar(@$aa) - 1)
	{
		my($point) = 0;
		while($point <= scalar(@aa_differences) - 1)
		{
			if($$aa[$k] eq $aa_differences[$point])
			{
				last;
			}else
			{
				$point += 1;
			}
		}
		if($point == scalar(@aa_differences))
		{
			$aa_differences[$pointer_to_dif] = $$aa[$k];
			$pointer_to_dif += 1;
		}
		$k += 1;
	}
	if(scalar(@aa_differences) > 1)
	{
		$k = 0;
		while($j <= scalar(@aa_differences) - 1)
		{
			$k = 0;
			$repetition[$j] = 0;
			while($k <= scalar(@$aa) - 1)
			{
				if($aa_differences[$j]  eq $$aa[$k])
				{
					$repetition[$j] += 1;
				}
				$k += 1;
			}
			$j += 1;
		}
		$k = 0;
		my($elements) = 0;
		while($k <= scalar(@repetition) - 1)
		{
		#	if($repetition[$k] == 1)
		#	{
		#		last;
		#	}else
		#	{
		#		$k += 1;
		#	}
			if($repetition[$k] > 1)
			{
				$elements += 1;
			}
			$k += 1;
		}
		#if($k == scalar(@repetition))
		if($elements > 1)
		{
			$$dif += 1;
		}
	}
}
sub Get_permutational_F
{
	my($correlations, $permut_prob) = @_;
	my($number) = 10000;
	for(my $i = 0; $i <= $number - 1; $i++)
	{
		$$permut_prob[$i] = $$correlations[int(rand(scalar(@$correlations)))];
	}
}
	