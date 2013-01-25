package InterSort;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;


#program by David to generate all co-evolving groups from CAPS Inter-molecular output
#created 21-3-2005
#altered 29-3-2005 to handle the new caps output file
#altered 9-6-2006 for inclusion as a sub-routien also added functionality to remove groups longer than a given threshold

sub make_final_groups 
{
	my ($first, $second) = @_;


	my @input = @$first;
	my $cutoff = $second;
	my $removing_groups = 0;
	#hash to store co-eveolver relationships
	my %data;


	foreach $_ (@input)
	{
		my $array_position = 0;
	
		my $site;
	
		my @current_data = split(/\s/, $_);
	
		foreach $_ (@current_data)
		{
		
			my $co_evolver = $_;
		
			if($array_position == 0)
			{
				$site = $_;
			}
		
			else
			{
				my $temp = $data{$co_evolver};
				$temp .= "$site,";
				$data{$co_evolver} = $temp;
				#print "$co_evolver -> $temp\n";
			}
		
			$array_position ++;
		}
	}

	my %new_hash;

	#go through the hash key by key
	foreach my $t1 (keys %data)           
	{                                     

		my $main_key = $t1;

		my $main_data = $data{$t1};

		my @site_holder;
	
		push (@site_holder, $main_key);


			#for each key go through the hash again and compare it to all other keys
			foreach my $t2 (keys %data)           
			{    
			
				#print OUTB "$main_key\t";
			
				my $current_key = $t2;
		
				my $current_data = $data{$t2};
	
				if (($current_key eq $main_key)||($current_data eq "done"))
				{
			
				}
		
				elsif ($current_data eq $main_data)
				{

					#print OUTB "$current_key\t";
					push (@site_holder, $current_key);
				
			
					$data{$main_key} = "done";
					$data{$current_key} = "done";
				}
	
			}

	
		############Switch Order#############
		my $primary_key;
		my $data;
		
		my @data_array = split(/,/, $main_data);
	
		foreach $_ (@data_array)
		{

			$primary_key .= "$_\t";
		}              

	
		foreach $_ (@site_holder)
		{
			$data .= "$_\t";
		}              
	
		$new_hash{$primary_key} = $data;
	}

 
	#calculate the total groups using pariwise comparisons
	#then compare the groups and remove redundant information
	#finally carry out the check that no group exceeds the specified length

	$removing_groups = 1;

	my %groups = make_groups(\%new_hash);
	my %reduced_groups = make_reduced_groups(\%groups);
	my %remove_list = get_group_lengths (\%reduced_groups, $cutoff);


	while ($removing_groups)
	{
		if ($remove_list{"passed"} == 1)
		{
			$removing_groups = 0;
		}

		else
		{
			%groups = remove_groups (\%groups, \%remove_list);
			%reduced_groups = make_reduced_groups(\%groups);
			%remove_list = get_group_lengths (\%reduced_groups, $cutoff);
		}
	}


	#returned hash contains sites from pritein 1 as the primary key pointing to their co-evolving sites from protein 2
	#(both are tab seperated strings of sites)
	return %reduced_groups
	
}


                                   
#takes a hash carries out pairwise comparison of groups and returns a hash of the newly created groups
sub make_groups
{
	my ($first) = @_;
        my %hash_in = %$first;
	my %hash_out;
	my $current_element = 0;
	foreach my $t1 (keys %hash_in)           
	{  

		my @sites = split(/\t/, $t1);
		my $temp = $hash_in{$t1};
		my @coevolvers = split(/\t/, $temp);
		my $length_sites;
	
	
		#go through each line in the file again to run the pairwise comparison 
		foreach my $t1 (keys %hash_in)           
		{  
			my @sites_B = split(/\t/, $t1);
			my $temp_B = $hash_in{$t1};
			my @coevolvers_B = split(/\t/, $temp_B);


			#run the comparison of the site arrays
			$length_sites = @sites;
			my $length_sites_B = @sites_B;
			my $number_matched = 0;
		
			#compare each element in the array
			foreach $_ (@sites)
			{
				$current_element = $_;	
				foreach $_ (@sites_B)
				{
					if ($current_element eq $_)
					{
						$number_matched++;
					}
				}
			}
		
			#if sites is a sub_set of sites_B
			if(($number_matched == $length_sites) && ($length_sites < $length_sites_B))
			{
				#add coevolvers_B to coevolvers
				foreach $_ (@coevolvers_B)
				{
				push (@coevolvers, $_);
				}
			}
		}
	
		$length_sites = @sites;
	
		if($length_sites == 0)
		{
		}
	
		else
		{
		
			my $catch_done = 0;
			my $key_string;
			my $data_string;
		
			#print the finalised output
			foreach $_ (@sites)
			{
				if ($_ eq "done")
				{
				$catch_done = 1;
				}
			
				elsif (!$catch_done)
				{
					$key_string .= "$_\t";
				}
			
			}
	
			if (!$catch_done)
			{
				foreach $_ (@coevolvers)
				{
					$data_string .= "$_\t";
				}
		
				$hash_out{$key_string} = $data_string;

			}
		
			$catch_done = 1;
		}
	}
return %hash_out;
}



#takes hash (as from make_reduced_groups) and an int and returns a list of entries that are too long
sub get_group_lengths
{
	my ($first, $second) = @_;
        my %hash_in = %$first;
	my %list_out;
	my $cutoff = $second;

        my $found_too_long = 0;
        
        $list_out{"passed"} = 0;
	#go through each key in the hash, if the number of elements in the key OR the data exceeds cutoff then add the key to the array
	
	foreach my $t1 (keys %hash_in)           
	{
		my $key = $t1;
	        my $data = $hash_in{$key};
	        my @key_list = split(/\t/, $key);
	        my @data_list = split(/\t/, $data);
        
	        my $key_length = @key_list;
	        my $data_length = @data_list;

	        if(($key_length > $cutoff) || ($data_length > $cutoff))
	        {
		        $found_too_long = 1;
		        $list_out{$key} = 1;
        	}
        
		else
        	{
        		$list_out{$key} = 0;
        	}

        }
        
        if (!$found_too_long)
        {
		$list_out{"passed"} = 1;
        }
        
	return %list_out;

}

#takes two hashes as input and removes from one those identified by the others
sub remove_groups
{

	my ($first, $second) = @_;
	my %hash_in = %$first;
	my %sites_to_remove = %$second;
	my %hash_out;

	foreach my $t1 (keys %hash_in)           
	{
		if ($sites_to_remove{$t1} == 0)
		{
			$hash_out{$t1} = $hash_in{$t1};
		}
	}

	return %hash_out;

}

#takes hash (as from make_groups) and through pairwise comparison removes overlapping data from groups
sub make_reduced_groups
{
#carry out pairwise comparisons

	my ($first) = @_;
	my %hash_in = %$first;
	my %sites_to_remove;
	my %hash_out;
	my $current_element = 0;

	#go through each line in the hash and store each key followed by its data followed by a zero
	foreach my $t1 (keys %hash_in)
	{ 
		my $key = $t1;
		my @data = split(/\t/, $hash_in{$key});
	
		foreach $_ (@data)
		{
		$sites_to_remove{$key}{$_} = 0;
		}
	}


	#go through each line in the hash
	foreach my $t1 (keys %hash_in)
	{  
		my @sites = split(/\t/, $t1);
		my $temp = $hash_in{$t1};
		my @coevolvers = split(/\t/, $temp);
		my $length_sites;
		my $key = $t1;
		my $is_subgroup = 0;

		#go through each line in the hash again to run the pairwise comparison
		foreach my $t1 (keys %hash_in)
		{
			my @sites_B = split(/\t/, $t1);
			my $temp_B = $hash_in{$t1};
			my @coevolvers_B = split(/\t/, $temp_B);
			my $key_B = $t1;



			#run the comparison of the site arrays
			$length_sites = @sites;
			my $length_sites_B = @sites_B;
			my $number_matched = 0;
		
			#compare each element in the array
			foreach $_ (@sites)
			{
				$current_element = $_;
				
				foreach $_ (@sites_B)
				{
					
					if ($current_element eq $_)
					{
						$number_matched++;
					}
				}
			
			}
			
			#if sites is a sub_set of sites_B
			if(($number_matched == $length_sites) && ($length_sites < $length_sites_B))
			{
				$is_subgroup = 1;

				#go through the coevolvers and mark the duplicates in the hash %sites_to_remove
				foreach $_ (@coevolvers)
				{
					my $element = $_;

					foreach $_ (@coevolvers_B)
					{
						my $element_B = $_;
						
						if($element eq $element_B)
						{
							$sites_to_remove{$key}{$element} = 1;
						}
					}

				}

			}
		}
	}

	#remove the duplicates from the input and send to output
	#go through each line in the hash
	foreach my $t1 (keys %hash_in)
	{  
		my $key = $t1;
		my $temp = $hash_in{$key};
		my @coevolvers = split(/\t/, $temp);

		my $replacement_data = "";

		#go through the coevolvers and choose which to remove
		foreach $_ (@coevolvers)
		{
			if($sites_to_remove{$key}{$_} == 0)
			{
				$replacement_data .= "$_\t"; 
			}
		}

		if ($replacement_data eq "")
		{}
		else
		{
			$hash_out{$key} = $replacement_data;
		}
	}
	
	return %hash_out;

}
1;