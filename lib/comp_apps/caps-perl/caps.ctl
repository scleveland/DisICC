Input file1: ~/Rails/MononegRails/DisorderConsensus/P_MAFFT_P_Cut_Fin_BEIV.fasta	* File containing sequence alignment for the first protein
Input file2: ~/Rails/MononegRails/DisorderConsensus/Para_L_Final_Padded_BEIV.fasta	* File containing sequence alignment for the second protein
Out file1: ~/Rails/MononegRails/DisorderConsensus/P_MAFFT_P_Cut_Fin_Para_L_Final_Padded_BEIV.out	* File where the output information should be stored
Co-evolution analysis: 1			* (0) Intra-molecular; (1) Inter-protein
Type of data 1: 0				* (0) amino acid alignment; (1) codon-based alignment
Type of data 2: 0				* (0) amino acid alignment; (1) codon-based alignment
3D test: 1						* Only applicable for intra-protein analysis: (0) perform test; (1) Test is not applicable
Reference sequence 1: 1	* the order in the alignment of the sequences giving the real positions in the protein for 3D analyses
Reference sequence 2: 1	* the order in the alignment of the sequences giving the real positions in the protein for 3D analyses
3D file: 				* name of the file containing 3D coordinates
Atom interval: 			* Amino acid atoms for which 3D coordinates are available
Significance test: 1				* (0) use threshold correlation; (1) random sampling
Threshold R-value: 0.1			* Threshold value for the correlation coeficient
Time correction: 1				* (0) no time correction; (1) weight correlations by the divergence time between sequences
Time estimation: 0				* (0) use synonymous distances by Li 1993; (1) use Poisson-corrected amino acid distances
Threshold alpha-value: 0.001		* This option valid only in case of random sampling
Random sampling: 10000			* Use in case significance test option is 1
Gaps: 2						* Remove all columns with gaps (0); Remove columns with a specified number of gaps (1); Do not remove columns with gaps(2)
Minimum R: 0.1					* Minimum value of correlation coeficient to be considered for filtering
GrSize: 3						* Maximum number of sites in the group permitted (given in percentage of protein length)
