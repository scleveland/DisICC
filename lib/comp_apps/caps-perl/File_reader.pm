package File_reader;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;
use Seq_manag;

sub Control_reader
{
	my($control, $information) = @_;
	my($Ctl) = Input_file($control);
	my($i) = 0;
	while(<$Ctl>)
	{
		$_ =~ s/ //g;
		my(@line) = split(/\:/, $_);
		my(@line2) = split(/\*/, $line[1]);
		$$information[$i] = $line2[0];
		$$information[$i] =~ s/\s+//g;
		$i += 1;
	}
}	 
sub Input_file
{
	my($infile) = @_;
	my($input_file);
	unless(open($input_file, $infile))
	{
		print STDERR "Cannot find $infile!!\n";
		exit;
	}
	return $input_file;
}

sub Output_file
{
	my($outfile) = @_;
	my($output_file);
	unless(open($output_file, ">$outfile"))
	{
		print STDERR "Cannot open $outfile!!\n";
		exit;
	}
	return $output_file;
}
sub Output_cvs
{
	my($outfile) = @_;
	my($output_file);
	unless(open($output_file, ">>$outfile"))
	{
		print STDERR "Cannot open $outfile!!\n";
		exit;
	}
	return $output_file;
}

sub Format_reader
{
	my($in, $Seqs, $names, $nuc_seq) = @_;
	my($infile) = Input_file($in);
	my($j, $longitud, $k, $controlador, $mega, $fasta, $phylip) = (0,0,0,0,0,0,0);
	my(@file_information);
	##Storing the file information into the array @file_information
	while(<$infile>)
	{
		chomp $_;
		$file_information[$j] = $_;
		$j += 1;
	}
	$j = 0;
	
	##Reading files in Phylip format
	if($file_information[0] =~ /^\d/)				##This will loop would read phylip format
	{
		while($j <= scalar(@file_information) - 1)
		{
			if(length($file_information[$j]) == 0)
			{
				$j += 1;
				next;
			}
			my(@array) = split(/\s+/, $file_information[$j]);
			if($j == 0)
			{
				$longitud = $array[1];
			}else
			{
				if($controlador == 0)
				{
				
					$$names[$k] = $array[0];
					$controlador += 1;
					if(scalar(@array) > 1)
					{
						$$Seqs[$k] = $array[1];
						$$nuc_seq[$k] = $array[1];
					}else
					{
						$$Seqs[$k] = '';
						$$nuc_seq[$k] = '';
					}
				}else
				{							#line 100
					if(length($$Seqs[$k]) == $longitud)
					{
						$k += 1;
						$$names[$k] = $array[0];
						if(scalar(@array) > 1)
						{
							$$Seqs[$k] = $array[1];
							$$nuc_seq[$k] = $array[1];
						}else
						{
							$$Seqs[$k] = '';
							$$nuc_seq[$k] = '';
						}
					}else
					{
						$$Seqs[$k] .= $array[0];
						$$nuc_seq[$k] = $array[0];
					}
				}
			}
			$j += 1;
		}
	}else ##Reading files in Mega or fasta format
	{
		while($j <= scalar(@file_information) - 1)
		{
			if(length($file_information[$j]) == 0)
			{
				$j += 1;
				next;
			}
			my(@array) = split(/\s+/, $file_information[$j]);
			if(($array[0] =~ /^(#mega)/i) || ($array[0] =~ /^(title:)/i))
			{
				$j += 1;
				next;
			}
			if($controlador == 0)
			{
				$$names[$k] = substr($array[0], 1, (length($array[0]) - 1));
				$controlador += 1;
				if(scalar(@array) > 1)
				{
					$$Seqs[$k] = $array[1];
					$$nuc_seq[$k] = $array[1];
				}else
				{
					$$Seqs[$k] = '';
					$$nuc_seq[$k] = '';
				}
			}else
			{
				if(($array[0] =~ /^>/) || ($array[0] =~ /^#/)) 
				{
					$k += 1;
					$$names[$k] = substr($array[0], 1, (length($array[0]) - 1));
					if(scalar(@array) > 1)
					{
						$$Seqs[$k] = $array[1];
						$$nuc_seq[$k] = $array[1];
					}else
					{
						$$Seqs[$k] = '';
						$$nuc_seq[$k] = '';
					}
				}else
				{
					$$Seqs[$k] .= $array[0];
					$$nuc_seq[$k] = $array[0];
				}
			}
			$j += 1;
		}
	}
	my($N_seqs) = 0;
	while($N_seqs <= scalar(@$Seqs) - 1)
	{
		chomp $$Seqs[$N_seqs];
		chomp $$nuc_seq[$N_seqs];
		$$Seqs[$N_seqs] =~ s/\r//g;
		$$nuc_seq[$N_seqs] =~ s/\r//g;
		$N_seqs += 1;
	}
}	
sub codon2aa
{
	my($sequences) = @_;
	my($i) = 0;
	while($i <= scalar(@$sequences) - 1)
	{
		my($seq) = $$sequences[$i];
		$$sequences[$i] = Seq_manag::nuc_amino($seq);
		$i += 1;
	}
}

sub test_seqlength
{
	my(@sequences) = @_;
	my($i) = 0;
	while($i <= scalar(@sequences) - 1)
	{
		if(length($sequences[$i]) != length($sequences[0]))
		{
			print "Error: Sequence " . ($i + 1) . " has different length\n";
			exit;
		}
		$i += 1;
	}
	return 0;
}
sub close_file
{
	my($file) = @_;
	close $$file;
}		
	
1;
