#include<iostream>
#include<fstream>
#include"BCFasta.h"
#include"file_manip.h"
#include"caps.h"
#include"create.h"
#include<getopt.h>
#include<math.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include<sys/time.h>
#include<iomanip>




#include <Seq/SequenceApplicationTools.h>
#include <Seq/SiteTools.h>
#include <Seq/SequenceTools.h>
#include <Seq/ioseq>
#include <Seq/alphabets>

// From Utils:
#include <Utils/ApplicationTools.h>
#include <Utils/TextTools.h>
#include <Utils/KeyvalTools.h>


// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/DRHomogeneousTreeLikelihood.h>
#include <Phyl/DRNonHomogeneousTreeLikelihood.h>
#include <Phyl/PatternTools.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/MarginalAncestralStateReconstruction.h>
#include <Phyl/OptimizationTools.h>
#include <Phyl/RASTools.h>
#include <Phyl/TreeLikelihoodTools.h>
#include <Phyl/DRTreeLikelihoodTools.h>
#include <Phyl/MarkovModulatedSubstitutionModel.h>
#include <Phyl/SubstitutionModelSet.h>
#include <Phyl/SubstitutionModelSetTools.h>
#include <Phyl/Newick.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/DataTable.h>
#include <NumCalc/MatrixTools.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/AutoParameter.h>

// From Utils:
#include <Utils/BppApplication.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>



//static int statl=0;

int main(int argc, char *argv[]){



	string pwd, temp_folder, mat_file, aminos, folder, forced, tree_file, mystring, structure, output, treefile, pdb_folder;
	double threshold, thresholdR, minR, gap_val, threshval=0.05;
	int num_randoms, three, time_corr, analysis=0, opt=0, time_type, seed, variable=0, tree_in=0, pdb_true=0;
	struct timeval start;
	const gsl_rng_type * T;
	gsl_rng *r;

	print_splash("stdout");

	if(argc<3){
		fprintf(stderr, "usage: [-F alignmentfolder] options([-S Structure folder]\n[--intra/--inter analysis type][-a alpha value][-r random samples][-N Newick formatted tree]\n");
		exit(-1);
	}

	Pwd(pwd);

	threshold=0.05;num_randoms = 100;thresholdR = 0.1;analysis=0;three=1;time_corr=0;

struct globalArgs_t {
    int noIndex;                /* -I option */
    char *langCode;             /* -l option */
    const char *outFileName;    /* -o option */
    FILE *outFile;
    int verbosity;              /* -v option */
    char **inputFiles;          /* input files */
    int numInputFiles;          /* # of input files */
    int randomized;             /* --randomize option */
} globalArgs;


static const struct option longOpts[] = {
    { "inter", no_argument, NULL, 'i' },
    { "intra", no_argument, NULL, 'l' },
    { NULL, no_argument, NULL, 0 }
};

int *longIndex;

	while ((opt = getopt_long(argc, argv, "H:b:c:F:o:r:ilg:a:G:N:vS:R:", longOpts, longIndex)) != -1) {
		switch(opt) {

			case 'N':
				tree_in = 1;
				treefile = optarg;
				break;
			case 'v':
				variable = 1;
				break;
			case 'i':
				analysis = 1;
				break;
			case 'l':
				analysis = 0;
				break;
			case 'a':
				threshval = atof(optarg);
				break;
			case 'S':
				pdb_true=1;
				pdb_folder = optarg;
				break;
			case 'F':
				if(optarg[0]=='~'||optarg[0]=='/'){
					mystring = optarg;
					break;
				}else{
					char temp[1000];
					getcwd(temp, 1000);
					mystring += temp;
					mystring  += "/";
					mystring += optarg;
					mystring += "/";
					break;
				}
			case 'R':
				thresholdR=atof(optarg);
				break;
			case 'r':
				num_randoms = atoi(optarg);	
				break;
			case 'g':
				gap_val=atof(optarg);
				break;
			case ':':
				fprintf(stderr, "Error: Unknown option passed in: %c\n", optopt); /* optarg defined in getopt.h */
		fprintf(stderr, "usage: [-F alignmentfolder] options([-S Structure folder]\n[--intra/--inter analysis type][-t threshold value][-r random samples][-N Newick formatted tree]\n");
				exit(1);
				break;
			case '?':
				fprintf(stderr, "Error: Unknown option passed in: %c\n", optopt);
		fprintf(stderr, "usage: [-F alignmentfolder] options([-S Structure folder]\n[--intra/--inter analysis type][-t threshold value][-r random samples][-N Newick formatted tree]\n");
				exit(1);
				break;
		}
	}

	Pwd(output);
	output += "/";

	if(mystring.length()==0){
		fprintf(stderr, "Error: you haven't entered an alignment folder! Use -F option.\n");
		exit(-1);
	}

	vector<string> files;
	files = Folder_to_vector(mystring.c_str());


	Fasta_vector file;
	Fasta_map file1, file2;

	/*set up random seed*/
	gettimeofday(&start, NULL);
	seed = start.tv_sec*start.tv_usec;
	gsl_rng_env_setup();
	T=gsl_rng_default;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, seed);
	/*finished initialising seed*/
	string B_aa="ACDEFGHIKLMNPQRSTVWY-";
	vector< vector<double> > D,D_corrected; 

	Blosum Blos62 = Blosum62();


	TreeTemplate<Node> *fixed = NULL;
	if(tree_in ==1){
		Newick * newickReader = new Newick(false);
		fixed = newickReader->read(treefile.c_str());
		delete newickReader;
	}

	vector<Node *> nams;


	ofstream interact;
	if(analysis==1){
		interact.open("coev_inter.csv");
		interact << "File1\tFile2\tCoevolving\tnum_pairs\ttotal_comp\tCut Off\tthreshold r\n";
	}

	vector<int> coev;
	vector<int> coev_total;
	double aver_data=0.0;

	vector< vector<double> > aver_back;
	vector<double> temp_back(files.size(), 0.0);
	for(int i=0;i<files.size();++i){
		aver_back.push_back(temp_back);
	}




	for(int i=0;i<files.size();++i){
		output.clear();
		Pwd(output);
		string temp_file;
		temp_file=mystring;
		temp_file+= "/";
		temp_file += files[i];
		vector<double> Correlations;




		if(analysis ==0){
			file.clear();
			Read_Fasta_vector(temp_file.c_str(), file);	
			file.check_for_pdb(pdb_folder, files[i]);
			file.check_valid_alignment("AMINO");


			vector<site> pairs_of_sites;
			output += files[i];
			output += ".out";
			cerr << "Input file: " << temp_file << endl;
			cerr << "Output will be sent to: " << output << endl << endl;
			print_splash(output);
			file.print_to_fasta(output.c_str());
			ofstream OUTPUT(output.c_str(), ios::app);
			int length = file.sequences[0].length();			

			/*get the JTT tree that corresponds to the sequences involved*/

			TreeTemplate<Node> *tree = NULL;

			std::auto_ptr<DistanceMatrix> DS;
			DS=ScoreDist(file.names, file.sequences);
			//			exit(-1);

			if(tree_in ==0){
				tree = create_input_tree(file.names, file.sequences);
			}else if(tree_in==1 && variable==1){

				string tempp = TreeTemplateTools::treeToParenthesis(*fixed, false);	
				tree = Remove_names(file.names, fixed);
				tempp = TreeTemplateTools::treeToParenthesis(*tree, false);
				vector<string> newnames = tree->getLeavesNames();
				nams = Put_distances_on_tree(tree, DS);


			}

			vector<double> distances;
			int k=1;
			int total=0;
			for(int j=0;j<file.names.size();++j){
				k=j+1;
				while(k<file.names.size()){
					//	cerr << "DS[" << j << "][" << k << "]=" << (*DS)(j, k) << endl;;
					distances.push_back((double)(*DS)(j,k));	
					//	cerr << "dist[" << total << "]=" << distances[total] << endl;
					++total;
					++k;
				}

			}

			vector<double> rel_dist = distances;
			sort(rel_dist.begin(), rel_dist.end());

			/*for(int j=0;j<distances.size();++j){
			  cerr << "dist[" << j << "]=" << distances[j] << endl;
			  distances[j]/=rel_dist[rel_dist.size()-1];
			  cerr << "dist[" << j << "]=" << distances[j] << endl;
			  }
			  */



			vector<double> totaltemp;
			double threshold;

			cerr << "Performing " << num_randoms << " simulations...\n\n";
			for(int j=0;j<num_randoms;++j){
				vector<string> Tnames, tsequences;
				create_seq(tree, Tnames, tsequences, length, "JC", variable, nams);
				//create_seq(tree, Tnames, tsequences, length, "JTT");

				vector< vector<double> > Dtemp, Dtempcorr;
				vector<int> temp_diff;
				Diff_aa_column(tsequences, temp_diff);
				Estimate_D(tsequences, Tnames, Dtemp, Dtempcorr, distances, time_corr, Blos62, temp_diff, tree);
				vector<double> tempCorrel;
				int numtempcor;

				numtempcor = intra(Dtemp, Dtempcorr, tempCorrel, time_corr,  (int)file.sequences[0].size(), 1);


				vector<double>::iterator it;
				it = tempCorrel.begin();
				totaltemp.insert(totaltemp.begin(), tempCorrel.begin(), tempCorrel.end());
				/*Do the same with JTT to get the number of pairs randomly for each alignment*/


			}

			/*Now do this with a better distribution*/
				cerr << "Performing Analysis...\n\n";

			for(int j=0;j<num_randoms;++j){
				vector<string> Tnames, tsequences;
				create_seq(tree, Tnames, tsequences, length, "JTT", variable, nams);
				//create_seq(tree, Tnames, tsequences, length, "JTT");

				vector< vector<double> > Dtemp, Dtempcorr;
				vector<int> temp_diff;
				Diff_aa_column(tsequences, temp_diff);
				Estimate_D(tsequences, Tnames, Dtemp, Dtempcorr, distances, time_corr, Blos62, temp_diff, tree);
				vector<double> tempCorrel;
				int numtempcor;

				numtempcor = intra(Dtemp, Dtempcorr, tempCorrel, time_corr,  (int)file.sequences[0].size(), 1);


				vector<double>::iterator it;
				it = tempCorrel.begin();
				totaltemp.insert(totaltemp.begin(), tempCorrel.begin(), tempCorrel.end());
				/*Do the same with JTT to get the number of pairs randomly for each alignment*/


				//cerr << "\r" << (float)j/(float)num_randoms << "\% finished";
				/*DO INTRA ON THESE*/
			}

			cerr << "\r100\% finished" << endl;
			//Get_threshold();
			sort(totaltemp.begin(), totaltemp.end());
			double meantemp = average_vec<double>(totaltemp);	
			int large = floor((totaltemp.size()*(1-threshval))) +1;
			threshold = totaltemp[large];

			/*get the number of columns with an average number of sites greater than the input value*/

			/*simulate a sequence alignment with the same distances as those given and a length equal to that of the input sequences*/
			/*run the intra coevolution analysis on this (random sample times) to get a cut off value*/

			vector< vector<double> > D, Dcorr;
			vector<int> aa_diff;
			Diff_aa_column(file.sequences, aa_diff);
			Estimate_D(file.sequences, file.names, D, Dcorr, distances, time_corr, Blos62, aa_diff, tree);
			vector<double> Correl;
			int numcor;
			numcor = intra(D, Dcorr, Correl, time_corr, (int)file.sequences[0].size(), 0);


				int ref=0;
				for(int p=0;p<file.names.size();++p){
				if(file.names[p]==file.ref_seq){
					ref=p;
					break;
				}	
				
				}

			print_to_file(Correl, D, Dcorr, threshold, file.sequences, aa_diff, OUTPUT, file, file.names[ref], files[i], ref, thresholdR);	
			OUTPUT.close();
			delete tree;
		}else{/*===============THIS IS THE INTER SECTION==============*/


			ofstream OUTPUT;
			file1.clear();



			Read_Fasta_map(temp_file.c_str(), file1);	
			/*loop over all files which are not the same as the i-th file*/
			for(int j=i+1;j<files.size();++j){
				output.clear();
				Pwd(output);
				file2.clear();
				Fasta_vector vec1, vec2;	
				string temp_file2;
				temp_file2 += mystring; 
				temp_file2 += "/";
				temp_file2 += files[j];
				Read_Fasta_map(temp_file2.c_str(), file2);
				Convert_to_vectors(file1, file2, vec1, vec2);	

				vec1.check_for_pdb(pdb_folder, files[i]);
				vec2.check_for_pdb(pdb_folder, files[j]);
				vec1.check_valid_alignment("AMINO");
				vec2.check_valid_alignment("AMINO");
				int sim_length1 = vec1.sequences[0].size() -vec1.getAveGaps();
				int sim_length2 = vec2.sequences[0].size() -vec2.getAveGaps();


				std::auto_ptr<DistanceMatrix> DS;
				vector<double> distances;
				//				if(i==0 && j==1){
				DS=ScoreDist(vec1.names, vec1.sequences);
				//exit(-1);
				int k=1;
				int total=0;
				for(int j=0;j<vec1.names.size();++j){
					k=j+1;
					while(k<vec1.names.size()){
						distances.push_back((double)(*DS)(j,k));	
						++total;
						++k;
					}

				}

				vector<double> rel_dist = distances;
				sort(rel_dist.begin(), rel_dist.end());
				//	}


				vector<site> pairs_of_sites;
				output += files[i];
				output += "_"; output += files[j];output += ".out";
				OUTPUT.open(output.c_str(), ios::app);
				cerr << "Input file1: " << files[i] << endl;
				OUTPUT << "Input file1: " << files[i] << endl;
				OUTPUT << "\n\nInput file2: " << files[j] << endl; 
				cerr << "\nInput file2: " << files[j] << endl; 
				cerr << "\nOutput will be sent to: " << output << endl << endl;
				interact << files[i] << "\t" << files[j] << "\t";


				print_splash(output);
				vec1.print_to_fasta(output.c_str());
				vec2.print_to_fasta(output.c_str());
				int length1 = vec1.sequences[0].length();			
				int length2 = vec2.sequences[0].length();

				TreeTemplate<Node> *tree1 = NULL;
				TreeTemplate<Node> *tree2 = NULL;


				if(tree_in ==0){
					tree1 = create_input_tree(vec1.names, vec1.sequences);
					tree2 = create_input_tree(vec2.names, vec2.sequences);
				}else if(tree_in ==1 && variable==1){

					string tempp = TreeTemplateTools::treeToParenthesis(*fixed, false);	
					tree1 = Remove_names(vec1.names, fixed);
					nams = Put_distances_on_tree(tree1, DS);

					tree2 = Remove_names(vec2.names, fixed);
					nams = Put_distances_on_tree(tree2, DS);

				}else{
					tree1 = fixed;
					tree2 = fixed;
				}	





				vector<double> totaltemp;
				double threshold;
				cerr << "Performing " << num_randoms << " simulations...\n\n";
				//				cerr << "num_randoms=" << num_randoms << endl;
				for(int j=0;j<num_randoms;++j){
					vector<string> Tnames1, tsequences1, Tnames2, tsequences2;
					create_seq(tree1, Tnames1, tsequences1, length1, "JC", variable, nams);
					create_seq(tree2, Tnames2, tsequences2, length2, "JC", variable, nams);

					vector< vector<double> > Dtemp1, Dtempcorr1, Dtemp2, Dtempcorr2;
					vector<int> temp_diff1, temp_diff2;
					Diff_aa_column(tsequences1, temp_diff1);
					Estimate_D(tsequences1, Tnames1,Dtemp1, Dtempcorr1, distances, time_corr, Blos62, temp_diff1, tree1);
					Diff_aa_column(tsequences2, temp_diff2);
					Estimate_D(tsequences2, Tnames2, Dtemp2, Dtempcorr2, distances, time_corr, Blos62, temp_diff2, tree1);
					vector<double> tempCorrel;
					inter(Dtemp1, Dtemp2, Dtempcorr1, Dtempcorr2, time_corr, tempCorrel, 1, temp_diff1, temp_diff2, tsequences1, tsequences2);

					vector<double>::iterator it;
					it = tempCorrel.begin();
					totaltemp.insert(totaltemp.begin(), tempCorrel.begin(), tempCorrel.end());
				}

				//tree1 = create_input_tree(vec1.names, vec1.sequences);
				vector<int> numbackground;
				sort(totaltemp.begin(), totaltemp.end());	

				int num_zero=0;
				for(int l=0;l<totaltemp.size();++l){
					if(totaltemp[l]==0.0){
						++num_zero;
						//	break;
					}
					//	cerr << "totaltemp[" << l << "]=" << totaltemp[l] << endl;
				}



				//	num_zero=0;
				//int value = floor(((totaltemp.size()-num_zero)*(1-(threshval/2))))+1+num_zero;
				int value = floor(((totaltemp.size())*(1-(threshval))))+1;


				//				cerr << "tottem_size=" << totaltemp.size() << endl;
				threshold = totaltemp[value];
				//exit(-1);
				/*double u = average_vec<double>(totaltemp);
				  double SDt = SD_vf(totaltemp, u);

				  cerr << "u=" << u << endl;
				  cerr << "SDt=" << SDt << endl;

				  threshold = gsl_cdf_gaussian_P(u, SDt);
				  */ 
				//				P_val.push_back(2*(1-gsl_cdf_gaussian_P((double)fabs(results[bottom]-mean), SD))); 
				//				cerr << "value=" << value << endl;
				//				cerr << "largest=" << totaltemp[0] << "\tsmallestl=" << totaltemp[totaltemp.size()-1] << endl; 
				//				cerr << "threshold=" << threshold << endl;
				/*Do the same with JTT to get the number of pairs randomly for each alignment*/




				int totjtt = 0;
				for(int j=0;j<num_randoms;++j){
					vector<string> Tnames1, tsequences1, Tnames2, tsequences2;
					Tnames1.clear();Tnames2.clear();tsequences1.clear();tsequences2.clear();
					vector< vector<double> > Dtemp1, Dtempcorr1, Dtemp2, Dtempcorr2;
					vector<int> temp_diff1, temp_diff2;
					create_seq(tree1, Tnames1, tsequences1, length1, "JTT", variable, nams);
					create_seq(tree2, Tnames2, tsequences2, length2, "JTT", variable, nams);
					Diff_aa_column(tsequences1, temp_diff1);
					Estimate_D(tsequences1, Tnames1, Dtemp1, Dtempcorr1, distances, time_corr, Blos62, temp_diff1, tree1);
					Diff_aa_column(tsequences2, temp_diff2);
					Estimate_D(tsequences2, Tnames2, Dtemp2, Dtempcorr2, distances, time_corr, Blos62, temp_diff2, tree1);
					vector<double> tempCorrel;
					/*DO INTER ON THESE*/
					inter(Dtemp1, Dtemp2, Dtempcorr1, Dtempcorr2, time_corr, tempCorrel, 0, temp_diff1, temp_diff2, tsequences1, tsequences2);
					totjtt += tempCorrel.size();
					numbackground.push_back(test_num(tempCorrel, threshold));
				}


				double average_back = average_vec<int>(numbackground);



				cerr << "Performing Analysis...\n\n";
				vector<double> Correl;
				vector< vector<double> > D, Dcorr, D2, Dcorr2;
				vector<int> diff1, diff2;
				Diff_aa_column(vec1.sequences, diff1);
				Diff_aa_column(vec2.sequences, diff2);
				Estimate_D(vec1.sequences, vec1.names, D, Dcorr, distances, time_corr, Blos62, diff1, tree1);
				Estimate_D(vec2.sequences, vec2.names, D2, Dcorr2, distances, time_corr, Blos62, diff2, tree1);
				inter(D, D2, Dcorr, Dcorr2, time_corr, Correl, 0, diff1, diff2, vec1.sequences, vec2.sequences);

				int ref1=0, ref2=0;
				for(int p=0;p<vec1.names.size();++p){
				if(vec1.names[p]==vec1.ref_seq){
					ref1=p;
					break;
				}	
				
				}

				
				for(int p=0;p<vec2.names.size();++p){
				if(vec2.names[p]==vec2.ref_seq){
					ref2=p;
					break;
				}	
				
				}



	
				int num_pairs = print_inter(Correl, threshold, OUTPUT, vec1.sequences, vec2.sequences, D, D2, diff1, diff2, vec1, vec2, vec1.names[ref1], vec2.names[ref2], output, files[i], files[j], ref1, thresholdR);

				double P_val=1.0;	
				int chival = Chi_squared(num_pairs, (int)Correl.size(), average_back, P_val);	
				if(chival==1){
					OUTPUT << "These proteins are possibly coevolving!!\tP_val<" << P_val << "\n\n";
					interact << "Yes\t" << num_pairs << "\t" << Correl.size() << "\t" << threshval << "\t" << threshold << endl;
				}else if(chival==2){
					OUTPUT << "These proteins have less than normal numbers of coevolving sites\tP_val<" << P_val << "\n\n";
					interact << "No-\t" << num_pairs << "\t" << Correl.size() << "\t" << threshval << "\t" << threshold << endl;
				}else{
					OUTPUT << "These proteins have a number of coevolvings sites which is similar to that expected for the distances involved\tP_val<" << P_val << "\n\n";
					interact << "No\t" << num_pairs << "\t" << Correl.size() << "\t" << threshval << "\t" << threshold << endl;
				}

				OUTPUT.close();

				//delete DS;
				delete tree2;
				//delete fixed;
				delete tree1;
				aver_data+=(double)num_pairs/(double)Correl.size();

				coev.push_back(num_pairs);
				coev_total.push_back((int)Correl.size());

				//	temp_aver.push_back((double)num_pairs/(double)Correl.size());
				aver_back[i][j] = (double)num_pairs/(double)Correl.size();
				aver_back[j][i] = (double)num_pairs/(double)Correl.size();
			}
		}		


		//	aver_back.push_back(temp_back);

	}

	interact << "\n\nInter data background normalisation of coevolution\n";
	interact << "=======================================================\n";

	if(analysis==1){
		aver_data/=(double)coev.size();


		map<string, double> back_map;

		int k=0;
		for(int i=0;i<files.size();++i){


			for(int j=i+1;j<files.size();++j){
				int nott = coev_total[k]-coev[k];
				double back_not = (1-aver_data)*(double)coev_total[k];
				double exp = aver_data*(double)coev_total[k];	
				double chi = pow(((double)coev[k]-exp),2)/exp;
				chi += pow(((double)nott-back_not),2)/back_not;

				interact << files[i] << "\t" << files[j];

				if(chi >=gsl_cdf_chisq_Pinv(0.95, 1) && (double)coev[k]>(exp)){
					interact << "\tYES\n";
				}else if(chi >=gsl_cdf_chisq_Pinv(0.95, 1) && (double)coev[k]<(exp)){
					interact << "\tNO-\n";
				}else{
					interact << "\tNO\n";
				}



				++k;
			}


		}


		interact << "\n\nWhat does each gene say about the others           \n";
		interact << "=======================================================\n";

		vector<double> file_back;
		for(int i=0;i<files.size();++i){

			double av =0.0;
			for(int j=0;j<files.size();++j){
				av += aver_back[i][j];
			}

			file_back.push_back(av/((double)files.size()-1));

		}

		k=0;
		for(int i=0;i<files.size();++i){

			for(int j=i+1;j<files.size();++j){


				int nott = coev_total[k]-coev[k];
				double back_not = (1-file_back[i])*(double)coev_total[k];
				double exp = file_back[i]*(double)coev_total[k];	
				double chi = pow(((double)coev[k]-exp),2)/exp;
				chi += pow(((double)nott-back_not),2)/back_not;

				interact << files[i] << "\t" << files[j];

				if(chi >=gsl_cdf_chisq_Pinv(0.95, 1) && (double)coev[k]>(exp)){
					interact << "\tYES\n";
				}else if(chi >=gsl_cdf_chisq_Pinv(0.95, 1) && (double)coev[k]<(exp)){
					interact << "\tNO-\n";
				}else{
					interact << "\tNO\n";
				}

				++k;

			}
		}


	}	


	interact.close();

	gsl_rng_free(r);
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Chi_squared
 *  Description:  Deduces the Probability of co-evolution of 2 proteins
 * =====================================================================================
 */
int Chi_squared (int num_pairs, int num_correlations, double background, double& P_val){

	double chi=0.0;

	//	chi +=	(pow((num_pairs-background),2))/background; 
	//	double nott=(num_correlations - num_pairs)/(double)num_correlations;
	//	double back_not = (num_correlations - background)/(double)num_correlations;
	double nott=(num_correlations - num_pairs);
	double back_not = (num_correlations - background);





	/* This is now the G-test, same theoretical assumptions */
	chi += (num_pairs*(log((double)num_pairs/(double)background)));
	chi += (nott*(log((double)nott/(double)back_not)));
	chi *= 2;

	//	cerr << "num_corr=" << num_correlations << "\tnum_pairs" << num_pairs << "\tbackgroud=" << background << endl;
	//	cerr << "chi=" << setprecision(14) << gsl_cdf_chisq_Pinv(gsl_cdf_chisq_P(chi, 1),1) << endl;

	if(chi>= gsl_cdf_chisq_Pinv(0.999, 1)){
		P_val = 0.001;
	}else if(chi>=gsl_cdf_chisq_Pinv(0.99, 1)){
		P_val = 0.01;
	}else if(chi>=gsl_cdf_chisq_Pinv(0.95, 1)){
		P_val = 0.05;
	}

	//	cerr << "chi=" << chi << endl;
	//	cerr << "PVAL=" << P_val << endl;

	if(chi >=gsl_cdf_chisq_Pinv(0.95, 1) && num_pairs>(background)){
		//	cerr << "linn=" << __LINE__ << endl;
		return 1;
	}else if(chi >=gsl_cdf_chisq_Pinv(0.95, 1) && num_pairs<(background)){
		//	cerr << "linn=" << __LINE__ << endl;
		return 2;
	}else{
		//	cerr << "linn=" << __LINE__ << endl;
		return 0;
	}


}		/* -----  end of function Chi_squared  ----- */





/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  print_inter
 *  Description:  Print the analysis results for inter molecular coevolution
 * =====================================================================================
 */
int print_inter(vector<double>& Correl, double threshold, ofstream& output, vector<string>& sequences, vector<string>& sequences2, vector< vector<double> >& D, vector< vector<double> >& D2, vector<int>& diff1, vector<int>& diff2, Fasta_vector& file1, Fasta_vector& file2, string& ref_spec1, string& ref_spec2, string& outfile, string& files1, string& files2, int ref, double& thresholdR){


	int gaps1=0;
	int cor=0, pairs=0;

	output << endl << endl;
	
	output << "Coevolving Pairs of amino acid sites\n";
	output << "=============================================================================\n";
	output << "Col1(real)\tCol2(real)\tDmean1\t\tDmean2\t\tCorrelation\n\n";
	output << "=============================================================================\n";

	double mean = average_vec<double>(Correl);
	double SD = SD_vf(Correl, mean);

	vector< vector<int> > pair;
	map< int, int> gap_map1, gap_map2;


	for(int i=0;i<(signed)sequences[0].size();++i){
		if(sequences[ref][i] == '-'){
			++gaps1;
		}
		int gaps2=0;
		for(int j=0;j<(signed)sequences2[0].size();++j){
			if(sequences2[ref][j]=='-'){
				++gaps2;	
			}

			if(fabs(Correl[cor])>=threshold && diff1[i]>0 && diff2[j]>0 && fabs(Correl[cor])>=thresholdR ){
				output << i+1 << "(" << i-gaps1+1 << ")\t\t" << j+1 << "(" << (j+1)-gaps2 << ")\t\t" << average_vec<double>(D[i]) << "\t\t" << average_vec<double>(D2[j]) << "\t\t" << Correl[cor] << endl; 
				++pairs;
				vector<int> tem;
				tem.push_back(i+1);
				tem.push_back(j+1);

				gap_map1[i+1]=i-gaps1+1;
				gap_map2[j+1]=j+1-gaps2;

				/*If structure info is present print distance between coevolving sites*/
				pair.push_back(tem);
			}	
			++cor;
		}
	}

	non_over_inter(pair, gap_map1, gap_map2, output, sequences, sequences2, file1, file2, ref_spec1, ref_spec2, outfile, files1, files2);
	group_inter(pair, gap_map1, gap_map2, output, sequences, sequences2, file1, file2, ref_spec1, ref_spec2, outfile, files1, files2);

	return pairs;
}		/* -----  end of function print_inter  ----- */



void non_over_inter (vector< vector<int> >& pairs, map<int, int >& gap_map1, map<int, int>& gap_map2, ofstream& output, vector<string>& sequences, vector<string>& sequences2, Fasta_vector& file1, Fasta_vector& file2, string& ref_spec,string& ref_spec2, string& outfile, string& files1, string& files2){



	
output << "Groups of coevolving amino acids" << endl;
output << "================================" << endl << endl;


	map<int, vector<int> > groups_left, groups_right;

	for(int i=0;i<pairs.size();++i){

		map<int, vector<int> >::iterator git;
		git=groups_left.find(pairs[i][0]);
		if(git==groups_left.end()){
			vector<int> tem;
			tem.push_back(pairs[i][1]);
			groups_left[pairs[i][0]]=tem;
		}else{
			git->second.push_back(pairs[i][1]);
		}

		git=groups_right.find(pairs[i][1]);
		
		if(git==groups_right.end()){
			vector<int> tem;
			tem.push_back(pairs[i][0]);
			groups_left[pairs[i][1]]=tem;
		}else{
			git->second.push_back(pairs[i][0]);
		}


	}



output << files1 << " vs " << files2 << "\n===========================================" << endl << endl;

ofstream fil1;
string tem; 
tem += files1;
tem += "_";
tem += files2;
tem += ".csv";

fil1.open(tem.c_str());

	fil1 << "sequenceName,msaColumnNumber,residueNumber,pdbNumber";

map<int, vector<int> >::iterator mit, mit2;
int p=1;
for(mit=groups_left.begin();mit!=groups_left.end();++mit){
	fil1 << ",Group" << p;
++p;	
}
fil1 << endl;
		


int k=1;
for(mit=groups_left.begin();mit!=groups_left.end();++mit){

	output << "Group " << k << ": " << mit->first << "(" << gap_map1[mit->first] << "):\t";

	for(int i=0;i<mit->second.size();++i){
		output << mit->second[i] << "(" << gap_map2[mit->second[i]] << ")\t"; 
		fil1 << ref_spec << "," << mit->second[i] << ",-1,-1";
			p=1;
			for(mit2=groups_left.begin();mit2!=groups_left.end();++mit2){

				if(p==k){
					fil1 << ",1";
				}else{
					fil1 << ",-1";
				}	
		
			}

			fil1 << endl;


	}
	++k;
	output << endl;
}

fil1.close();


output << "File 2 vs File 1" << "===========================================" << endl << endl;

tem.clear();
tem += files2;
tem += "_";
tem += files1;
tem += ".csv";


k=1;
for(mit=groups_right.begin();mit!=groups_right.end();++mit){


	output << "Group " << k << ": " << mit->first << "(" << gap_map2[mit->first] << "):\t";

	for(int i=0;i<mit->second.size();++i){
		output << mit->second[i] << "(" << gap_map1[mit->second[i]] << ")\t"; 

		fil1 << ref_spec << "," << mit->second[i] << ",-1,-1";
			p=1;
			for(mit2=groups_left.begin();mit2!=groups_left.end();++mit2){

				if(p==k){
					fil1 << ",1";
				}else{
					fil1 << ",-1";
				}	
		
			}

			fil1 << endl;



	}
	output << endl;
	++k;
}
output << endl;




}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  group_inter
 *  Description:  
 * =====================================================================================
 */
void group_inter (vector< vector<int> >& pairs, map<int, int >& gap_map1, map<int, int>& gap_map2, ofstream& output, vector<string>& sequences, vector<string>& sequences2, Fasta_vector& file1, Fasta_vector& file2, string& ref_spec,string& ref_spec2, string& outfile, string& files1, string& files2){


	map<int, vector<int> > groups;
	map<int, vector<int> > groups2;
	for(int i=0;i<pairs.size();++i){

		map<int, vector<int> >::iterator it;
		it = groups.find(pairs[i][0]);
		if(it==groups.end()){
			vector<int> tem;
			tem.push_back(pairs[i][1]);
			groups[pairs[i][0]]=tem;
		}else{
			it->second.push_back(pairs[i][1]);
		}
	}


	for(int i=0;i<pairs.size();++i){

		map<int, vector<int> >::iterator it;
		it = groups2.find(pairs[i][1]);
		if(it==groups2.end()){
			vector<int> tem;
			tem.push_back(pairs[i][0]);
			groups2[pairs[i][1]]=tem;
		}else{
			it->second.push_back(pairs[i][0]);
		}
	}

	output << endl << "Overlapping groups of coevolving residues" << endl << "====================================================================" << endl;
	output << "Group #: Coevolving Sites" << endl << endl;

	map<int, vector<int> >::iterator git, git2;


	map< vector<int>, vector<int> > full_map;
	map< vector<int>, vector<int> >::iterator mit;;




	vector< vector<int> > left, right;

	for(git=groups.begin();git!=groups.end();++git){
		vector<int> temp;
		temp.push_back(git->first);
		left.push_back(temp);
		right.push_back(git->second);
	}


	for(git=groups2.begin();git!=groups2.end();++git){
		for(int i=0;i<left.size();++i){
			vector<int>::iterator vit, vit2;

			vit = find(right[i].begin(), right[i].end(), git->first);
			if(vit!=right[i].end()){
				for(int j=0;j<git->second.size();++j){
					vit2 = find(left[i].begin(), left[i].end(), git->second[j]);
					if(vit2==left[i].end()){
						vector<int> v(10000);
						vector<int>::iterator it;
						it = set_union (left[i].begin(), left[i].end(), git->second.begin(), git->second.end(), v.begin());
						//			cerr << "numsize=" << int(it - v.begin()) << endl;
						v.erase(v.begin()+int(it -v.begin()), v.end());

						//	left[i].push_back(git->second[j]);
						left[i].clear();
						left[i]=v;
					}
				}
			}


		}
	}


	for(git=groups.begin();git!=groups.end();++git){
		for(int i=0;i<left.size();++i){
			vector<int>::iterator vit, vit2;

			vit = find(left[i].begin(), left[i].end(), git->first);
			if(vit!=left[i].end()){
				for(int j=0;j<git->second.size();++j){
					vit2 = find(right[i].begin(), right[i].end(), git->second[j]);
					if(vit2==right[i].end()){
						vector<int> v(10000);
						vector<int>::iterator it;
						it = set_union (right[i].begin(), right[i].end(), git->second.begin(), git->second.end(), v.begin());

						v.erase(v.begin()+int(it -v.begin()), v.end());
						right[i].clear();
						right[i]=v;

					}

				}
			}


		}
	}


	int k=1;
	map< vector<int>, vector<int> > mymap;
	map< vector<int>, vector<int> >::iterator myit, myit2;


	for(int i=0;i<right.size();++i){
		sort(left[i].begin(), left[i].end());
		sort(right[i].begin(), right[i].end());
		mymap[left[i]]=right[i];
	}

	int brekkie=1;

	for(myit=mymap.begin();myit!=mymap.end();++myit){
		if(brekkie==1){
			myit=mymap.begin();
			brekkie=0;
		}


		myit2=myit;
		++myit2;
		for(;myit2!=mymap.end();++myit2){
			vector<int> v(10000), v2(10000);
			vector<int>::iterator it, it2;
			it=set_intersection (myit->first.begin(), myit->first.end(), myit2->first.begin(), myit2->first.end(), v.begin());
			it2=set_intersection (myit->second.begin(), myit->second.end(), myit2->second.begin(), myit2->second.end(), v2.begin());

			if(int(it -v.begin())>0 || int(it2 - v2.begin())>0){
				vector<int> templ(20000), tempr(20000);
				v.erase(v.begin()+int(it -v.begin()), v.end());
				v2.erase(v2.begin()+int(it -v2.begin()), v2.end());

				it=set_union(myit->first.begin(), myit->first.end(), myit2->first.begin(), myit2->first.end(), templ.begin());
				it2=set_union(myit->second.begin(), myit->second.end(), myit2->second.begin(), myit2->second.end(), tempr.begin());
				mymap.erase(myit);
				templ.erase(templ.begin()+int(it -templ.begin()), templ.end());
				tempr.erase(tempr.begin()+int(it2 -tempr.begin()), tempr.end());
				mymap.erase(myit2);
				mymap.insert (mymap.begin(), pair<vector<int> ,vector<int> >(templ,tempr));
				myit=mymap.begin();
				brekkie=1;
				break;	
			}
		}

	}



	ofstream gro, gro2;
	string tem;
	tem += files1;
	tem += "_";
	tem += files2;
	tem += "-overlap.csv";
cerr << "Residue annotations sent to: " << tem << endl << " and : ";
	gro.open(tem.c_str());
	tem.clear();
	tem += files2;
	tem += "_";
	tem += files1;
	tem += "-overlap.csv";
cerr << tem << endl;
	gro2.open(tem.c_str());



	gro << "sequenceName,msaColumnNumber,residueNumber,pdbNumber,score,Group\n";
	gro2 << "sequenceName,msaColumnNumber,residueNumber,pdbNumber,score,Group\n";

	for(myit=mymap.begin();myit!=mymap.end();++myit){

		output << "Group " << k << ": [";
		for(int i=0;i<myit->first.size();++i){
			output << myit->first[i] << "(" << gap_map1[myit->first[i]] << ")\t";
		}


		output << " | " ;
		for(int i=0;i<myit->second.size();++i){
			output << myit->second[i] << "(" << gap_map2[myit->second[i]] << ")\t";
		}

		output << " ] \n";
		if(file1.has_pdb==true && file2.has_pdb==true){
			output << "[ ";
			vector<int> temp = myit->first;
			output << group_dist_inter(temp, 0, file1, gap_map1, ref_spec, gro, k, temp.size());
			output << " | ";

			output << group_dist_inter(myit->second, 0, file2, gap_map2, ref_spec2, gro2, k, myit->second.size());
			output << " ]\n";

		}

		output << "==========================================================\n";

		++k;
	}

gro.close();
gro2.close();



	/*	int k=1;
		for(git=groups.begin();git!=groups.end();++git){
		output << "Group " << k << ": " << git->first << " [ ";
		vector<int>::iterator vit;
		for(vit=git->second.begin();vit!=git->second.end();++vit){
		output << *vit << "\t";
		}
		output << " ]" << endl;
		++k;

		}	
		*/






}		/* -----  end of function group_inter  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Network::add_connection
 *  Description:  
 * =====================================================================================
 */
void Network::add_connection (int& i, int& j){


/*find i, if not there add it*/
		map<int, int>::iterator mit1;
		mit1 = node_map.find(i);

	if(mit1==node_map.end()){
		


	Nod temp;
	vector<int> tem1, tem2;
	tem1.push_back(j);
	temp.node_connect[i]=tem1;
	nodes.push_back(temp);
	node_map[i]=node_map.size()-1;;

	}else{
		
		vector<int>::iterator vit;
	
int tes = mit1->second;
//vector<int> tem_vec = nodes[tes].node_connect->second;
		//ves = nodes[tes].node_connect->second();	
		vit = find(nodes[tes].node_connect[tes].begin(), nodes[tes].node_connect[tes].end(), j);
		//map<int, int>::iterator mit;
		//mit = node_map.find(j);

		if(vit!=nodes[tes].node_connect[tes].end()){
			nodes[tes].node_connect[tes].push_back(j);
		}

	}

mit1 = node_map.find(j);
	if(mit1==node_map.end()){
	Nod temp;
	vector<int> tem1, tem2;
	tem2.push_back(i);
	temp.node_connect[i]=tem1;
	nodes.push_back(temp);
	node_map[j]=node_map.size()-1;

	}else{
		vector<int>::iterator vit;
int tes = mit1->second;
		vit = find(nodes[mit1->second].node_connect[tes].begin(), nodes[mit1->second].node_connect[tes].end(), i);
		if(vit!=nodes[mit1->second].node_connect[tes].end()){
			nodes[mit1->second].node_connect[tes].push_back(i);
		}
	}

}		/* -----  end of function Network::add_connection  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Network:print_network
 *  Description:  
 * =====================================================================================
 */
void Network::print_network (ofstream& file){

for(int i=0;i<nodes.size();++i){
	file << "Group" << i << " :";

	map<int, vector<int> >::iterator it;

	for(it=nodes[i].node_connect.begin();it!=nodes[i].node_connect.end();++it){
		file << it->first << "\t";
		for(int j=0;j<it->second.size();++j){
			file << it->second[j] << "\t";
		}
	file << endl;

	}

}

}		/* -----  end of function Network:print_network  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  test_num
 *  Description:  test how many values are greater than the threshold 
 * =====================================================================================
 */
int test_num (vector<double>& Correl, double thresh){

	int number=0;

	for(int i=Correl.size()-1;i>=0;--i){
		//cerr << "Correl[" << i << "]=" << Correl[i] << "thresh=" << thresh << endl;
		if(fabs(Correl[i])>=fabs(thresh)){
			++number;
		}
	}

	return number;
}		/* -----  end of function test_num  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  inter
 *  Description:  Inter Molecular Co-evolution analysis
 * =====================================================================================
 */
void inter(vector< vector<double> >& D, vector< vector<double> >& D2, vector< vector<double> >& D_correct, vector< vector<double> >& D_correct2, int& correction, vector<double>& Correl, int simulate, vector<int>& diff1, vector<int>& diff2, vector<string>& sequences, vector<string>& sequences2){


	for(int i=0;i<D.size();++i){
		vector<double> D_column(D[0].size(), 0.0), D_column2(D[0].size(), 0.0);
		//		D_column.assign(D[0].size(), 0.0);
		//		D_column2.assign(D[0].size(), 0.0);
		if(correction == 0)
		{
			for(int m=0;m<D[0].size();++m){
				D_column[m] = D[i][m];
			}
		}else
		{
			for(int m=0;m<D[0].size();++m){
				D_column[m] = D_correct[i][m];
			}
		}
		for(int j=0;j<D2.size();++j){
			if(diff2[j]==1 && diff1[i]==1){
				if(correction == 0)
				{
					for(int n=0;n<D[0].size();++n){
						D_column2[n] = D2[j][n];
					}
				}else
				{
					for(int n=0;n<D[0].size();++n){
						D_column2[n] = D_correct2[j][n];
					}
				}

				if(simulate==1){
					//Correl.push_back(fabs(Correlation(D_column, D_column2)));
					/*if(Correlation(D_column, D_column2)==1){
					//cerr << "atat=" << atanh(0.999999);
					Correl.push_back(atanh(0.9999999999));
					}else{	*/
					Correl.push_back(fabs(Correlation(D_column, D_column2)));
					//	}
				}else{
					//						cerr << "atanh=" << atanh(Correlation(D_column, D_column2)) << "\tCorrel=" << Correlation(D_column, D_column2) << endl;
					/*	if(Correlation(D_column, D_column2)==1){
						cerr << "atat=" << atanh(0.99999999999) << endl;;
						Correl.push_back(atanh(0.9999999999));
						}else{*/	
					//if(check_gaps(sequences, sequences2, i, j)){					
					Correl.push_back(Correlation(D_column, D_column2));
					//	}else{
					//		Correl.push_back(0.0);
					//	}
					//		}
				}
			}else{
				//					cerr << "atanh=" << atanh(Correlation(D_column, D_column2)) << "\tCorrel=" << Correlation(D_column, D_column2) << endl;
				//	Correl.push_back(Correlation(D_column, D_column2));
				Correl.push_back(0.0);
				//	if(simulate==0)
				//	cerr << "ii=" << i << "\tj=" << j << endl;
			}
		}/*end while j*/
	}/*end while i*/
}		/* -----  end of function inter  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  check_gaps
 *  Description:  
 * =====================================================================================
 */
bool check_gaps (vector<string>& seq, vector<string>& seq2, int i, int j)
{

	bool test=true;
	int count=0;
	for(int k=0;k<seq.size();++k){
		if(seq[k][i]=='-'){
			count++;
		}
	}

	if((double)count<(0.8*seq.size())){
		test=false;
	}

	if(count>=10){
		test=true;
	}

	count =0;
	for(int k=0;k<seq2.size();++k){
		if(seq2[k][i]=='-'){
			count++;
		}
	}

	if((double)count<(0.8*seq.size())){
		test=false;
	}

	if(count>=10){
		test=true;
	}

	return test;;
}		/* -----  end of function check_gaps  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Convert_to_vectors
 *  Description:  For inter-coev analysis, converts maps to vectors of common species 
 * =====================================================================================
 */
void Convert_to_vectors(Fasta_map& file1, Fasta_map& file2, Fasta_vector& vec1, Fasta_vector& vec2){

	map<string, string>::iterator it, fin;

	vec1.filename = file1.filename;
	vec2.filename = file2.filename;
vec1.has_pdb=false;
vec2.has_pdb=false;
	for(it=file1.sequences.begin();it!=file1.sequences.end();++it){

		fin = file2.sequences.find(it->first);
		if(fin!=file2.sequences.end()){
			vec1.names.push_back(it->first);
			vec1.sequences.push_back(it->second);
			vec2.names.push_back(fin->first);
			vec2.sequences.push_back(fin->second);
		}
	}
}		/* -----  end of function Convert_to_vectors  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  print_to_file
 *  Description:  looks for the Correlations, tests if they are greater than threshold
 * =====================================================================================
 */

void print_to_file (vector<double>& Correl, vector< vector<double> >& D, vector< vector<double> >& D_corr, double threshold, vector<string>& sequences, vector<int>& aa_differences, ofstream& output, Fasta_vector& file, string& ref_spec, string& filesi, int ref, double& thresholdR){
	int gaps1=0, cor=0;




	output << "\nCoevolution Analysis\n";
	output << "\n=======================================================================\n";
	output <<  "Position AA1\tPosition AA2\tMean D1\t\tMean D2\t\tCorrelation\n";
	output <<  "------------\t------------\t-------\t\t-------\t\t-----------\n";

	map<int, int> gap_map;

	vector< vector<int> > pairs_of_sites;
	//cerr << "seq_length=" << sequences[0].size() << endl;
	for(int i=0;i<(signed)sequences[0].size()-1;++i){
		if(sequences[ref][i]=='-'){
			++gaps1;
		}

		int gaps2=0;	
		for(int j=i+1;j<(signed)sequences[0].size();++j){
			if(sequences[ref][j]=='-'){
				++gaps2;
			}	

			if(fabs(Correl[cor])>=threshold && aa_differences[i]>0 && aa_differences[j]>0 && fabs(Correl[cor])>=thresholdR){
				output << i+1 << "(" << i+1-gaps1 << ")\t\t" << j+1 << "(" << j+1-(gaps1+gaps2) << ")\t\t" << average_vec<double>(D[i]) << "\t\t" << average_vec<double>(D[j]) << "\t\t" << Correl[cor];
				vector<int> tem;
				tem.push_back(i+1);
				tem.push_back(j+1);
				output << "\t" << assess_hyd(sequences, i, j) << "\t";
				gap_map[i+1]=i+1-gaps1;
				gap_map[j+1]=j+1-(gaps1+gaps2);

				if(file.has_pdb==true){
					int real_pos1=i-gaps1;
					int real_pos2=j-(gaps1+gaps2);
					output << get_distance(real_pos1, real_pos2, file) << "\t";
				}

				output << endl;

				pairs_of_sites.push_back(tem);
			}

			++cor;

		}


	}


	non_over_intra(pairs_of_sites, output, file, gap_map, ref_spec, filesi);
	group_intra(pairs_of_sites, output, file, gap_map, ref_spec, filesi);
	//	Hydro(pairs_of_sites, sequences, output);

}		/* -----  end of function print_to_file  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  get_distance
 *  Description:  
 * =====================================================================================
 */
double  get_distance (int pos1, int pos2, Fasta_vector& file)
{

	double dist;


	map<int, int>::iterator mit;

	/*for(mit=file.pdb.chains[0].res_map.begin();mit!=file.pdb.chains[0].res_map.end();++mit){

	  cerr << "pos=" << mit->first << "\tres_map=" << mit->second << endl;
	  }
	  */
	//exit(-1);

	//cerr << "res_map[" << pos1 << "]=" << file.pdb.chains[0].res_map[pos1] << "\tres_map[" << pos2 << "]=" << file.pdb.chains[0].res_map[pos2] << endl;
	//		double x = pow(file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos1]].mean_pos.x - file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos2]].mean_pos.x, 2);
	//		double y = pow(file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos1]].mean_pos.y -file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos2]].mean_pos.y, 2);
	//		double z = pow(file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos1]].mean_pos.z -file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos2]].mean_pos.z, 2);
	//	cerr << "x=" << x << "\ty=" << y << "\tz=" << z << endl;
	if(file.pdb.chains[0].res_map.find(pos1)==file.pdb.chains[0].res_map.end() || file.pdb.chains[0].res_map.find(pos2)==file.pdb.chains[0].res_map.end()){
		return 999.9;
	}else{
		dist = sqrt((pow(file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos1]].mean_pos.x - file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos2]].mean_pos.x, 2)) +(pow(file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos1]].mean_pos.y -file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos2]].mean_pos.y, 2) )+(pow(file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos1]].mean_pos.z -file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[pos2]].mean_pos.z, 2) ) );

	}
	return dist;
}		/* -----  end of function get_distance  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  assess_hyd
 *  Description:  
 * =====================================================================================
 */
string assess_hyd(vector<string>& sequences, int col1, int col2){

	hydro myhyd;
	myhyd.set_hyd();
	string str;


	double meanh=0.0;
	/*first get average hydro and MW*/
	for(int i=0;i<sequences.size();++i){
		meanh+= myhyd.values[sequences[i][col1]];					
		meanh+= myhyd.values[sequences[i][col2]];					
	}

	meanh /=(sequences.size()*2);

	double max = (155*0.05);
	double lower = meanh - max;
	double upper = meanh + max;
	int t1=0;
	for(int i=0;i<sequences.size();++i){
		double tempd=myhyd.values[sequences[i][col1]];
		tempd+=myhyd.values[sequences[i][col1]];
		tempd/=2;
		if(tempd<lower || tempd>upper){
			t1=1;
		}

	}
	MW mymw;
	mymw.set_mw();

	meanh=0.0;
	for(int i=0;i<sequences.size();++i){
		meanh+= mymw.values[sequences[i][col1]];					
		meanh+= mymw.values[sequences[i][col2]];					
	}

	meanh /=(sequences.size()*2);

	max = (129*0.05);
	lower = meanh - max;
	upper = meanh + max;
	int t2=0;


	for(int i=0;i<sequences.size();++i){
		double tempd=mymw.values[sequences[i][col1]];
		tempd+=mymw.values[sequences[i][col1]];
		tempd/=2;
		if(tempd<lower || tempd>upper){
			t2=1;
		}
	}


	if(t1==0 && t2==0){
		str += "HYD & MW";
	}else if(t1==0){
		str += "HYD";
	}else if(t2==0){
		str += "MW";
	}else{
		str += "NO INFO";
	}
	return str;
}		/* -----  end of function assess_hyd  ----- */





/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Hydro
 *  Description:  
 * =====================================================================================
 */
void Hydro (vector< vector<int> >& pairs, vector<string>& sequences, ofstream& output)
{
	hydro myhyd;
	myhyd.set_hyd();


	for(int i=0;i<pairs.size();++i){
		//		get_average_hyd		




	}




}		/* -----  end of function Hydro  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  non_over_intra
 *  Description:  
 * =====================================================================================
 */
	void non_over_intra (vector< vector<int> >& pairs, ofstream& output, Fasta_vector& file, map<int, int>& gap_map, string& ref_spec, string& filei){

	map<int, vector<int> > map_group;

	int last = pairs[0][0];
	vector<int> tem_vec;
	tem_vec.push_back(pairs[0][1]);
	map_group[pairs[0][0]] = tem_vec;


/* make a map of each residue to its coevolving sites*/
	for(int i=1;i<pairs.size();++i){

		if(pairs[i][0]==last){
			map_group[pairs[i][0]].push_back(pairs[i][1]);
		}else{
			last = pairs[i][0];
			tem_vec.clear();
			tem_vec.push_back(pairs[i][1]);
			map_group[last]=tem_vec;
		}

	}



	map<int, vector<int> >::iterator mit;

ofstream fil1;
string tem;
tem += filei;
tem += ".csv";
fil1.open(tem.c_str());

vector< vector<int> > total_group;


output << endl << "Groups of coevolving Amino Acids" << endl << "==========================================\n";

	int k=1;

vector< vector<int> > subs;

	for(mit=map_group.begin();mit!=map_group.end();++mit){
//	cerr << mit->first << ":\t";
		for(int i=0;i<mit->second.size();++i){
	//		cerr << mit->second[i] << "\t";
			vector<int> tem;
				tem.push_back(mit->first);
				tem.push_back(mit->second[i]);
			for(int j=i+1;j<mit->second.size();++j){
				map<int, vector<int> >::iterator git, git2;
				git = map_group.find(mit->second[i]);
			//	git2 = map_group.find(mit->second[j]);
			
			if(git!=map_group.end()){
					vector<int>::iterator vit;
					//vit = git->second.find(mit->second[j]);
					vit = find(git->second.begin(), git->second.end(), mit->second[j]);
					if(vit!=git->second.end()){
						tem.push_back(mit->second[j]);
					}
			}	
				
			}
			
			if(tem.size()>2 && not_subset(subs, tem)) {
				output << "Group " << k << ":\t";
				for(int l=0;l<tem.size();++l){

				output << tem[l]  << "(" << gap_map[tem[l]] << ")\t"; 
//					gro << ref_spec << "," << sites[i] << ",-1,-1," << "," << k << endl; 
				}
			total_group.push_back(tem);
				output << endl;
output << "=========================================================\n";
				++k;
				subs.push_back(tem);
			}



		}
//	cerr << endl;
	}
output << endl;

cerr << "Residue annotations sent to: " << tem << endl;
	
	fil1 << "sequenceName,msaColumnNumber,residueNumber,pdbNumber";
for(int i=0;i<total_group.size();++i){
	fil1 << ",Group" << i+1;
}	
fil1 << endl;



for(int i=0;i<total_group.size();++i){
	for(int j=0;j<total_group[i].size();++j){
					fil1 << ref_spec << "," << total_group[i][j] << ",-1,-1";
					for(int l=0;l<total_group.size();++l){
						if(l==i){
							fil1 << ",1";
						}else{
							fil1 << ",-1";
						}
					}
					fil1 << endl;
 
	}
}



	fil1.close();
	



}		/* -----  end of function non_over_intra  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  is_subset
 *  Description:  
 * =====================================================================================
 */
	bool not_subset (vector< vector<int> >& subs, vector<int>& set ){



		for(int i=0;i<subs.size();++i){
	int hits=0;
			
			for(int j=0;j<set.size();++j){
				if(find(subs[i].begin(), subs[i].end(), set[j])!=subs[i].end()){
					++hits;
				}
			}

		if(hits==set.size()){
			return false;
		}


		}

	

	return true;
}		/* -----  end of function is_subset  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  group_intra
 *  Description:  
 * =====================================================================================
 */
void group_intra (vector< vector<int> >& pairs, ofstream& output, Fasta_vector& file, map<int, int>& gap_map, string& ref_spec, string& filei){


int cur = pairs[0][0];
vector< vector<int> > groups;
vector<int> tg;
tg.push_back(pairs[0][0]);

	for(int i=0;i<pairs.size();++i){
		if(pairs[i][0]==cur){
			tg.push_back(pairs[i][1]);		
		}else{
		cur = pairs[i][0];
		groups.push_back(tg);
		tg.clear();
		tg.push_back(pairs[i][0]);
		tg.push_back(pairs[i][1]);
	
		}
	

	}


int reset=0;
for(int i=0;i<groups.size();++i){


	for(int j=i+1;j<groups.size();++j){
		
		if(reset==1){
		i=0;
		j=1;
		reset=0;
		}

		vector<int> mer(groups[i].size()+groups[j].size());
		merge(groups[i].begin(), groups[i].end(),groups[j].begin(),groups[j].end(), mer.begin());

			map<int, int> mymap;
			for(int k=0;k<mer.size();++k){
				mymap[mer[k]]=1;
			} 
			if(mymap.size()<(groups[i].size()+groups[j].size())){
				groups.erase(groups.begin()+i);
				groups.erase(groups.begin()+j);
				vector<int> tem;
				map<int, int>::iterator mit;

				for(mit=mymap.begin();mit!=mymap.end();++mit){
					tem.push_back(mit->first);
				}
				groups.push_back(tem);
				reset=1;

			}
	}
}




output << endl << endl << "Overlapping groups of coevolving amino acids" << endl << "=================================================" << endl;

	ofstream gro;
	string tem;
	tem += filei;
	tem += "-overlap.csv";
cerr << "Over-lapping residue annotations sent to: " << tem << endl;

	gro.open(tem.c_str());
	gro << "sequenceName,msaColumnNumber,residueNumber,pdbNumber,score,Group\n";

for(int i=0;i<groups.size();++i){
	output << "Group " << i+1 << ": ";
	for(int j=0;j<groups[i].size();++j){
		output << groups[i][j] << "(" << gap_map[groups[i][j]] << ")\t"; 
	}

		if(file.has_pdb==true){
			output << "mean distance: " << group_dist_intra(groups[i], 0, file, gap_map, ref_spec, gro, i, groups[i].size());
		}


	output << endl;
	output << "==================================================" << endl;
}

gro.close();



}		/* -----  end of function group_intra  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  group_dist
 *  Description:  
 * =====================================================================================
 */
double group_dist_intra (vector<int>& sites, int site, Fasta_vector& file , map<int, int>& gap_map, string& species, ofstream& gro, int& group_num, unsigned int group_size)
{
	double mean=0.0;
	int k=0;

	for(int i=0;i<sites.size();++i){
		for(int j=i+1;j<sites.size();++j){	
			if(check_existance(gap_map[sites[i]], gap_map[sites[j]], file)){
				double temp = get_distance(gap_map[sites[i]], gap_map[sites[j]], file);
				
				if(temp!=999.9){
					mean += temp;
					++k;
//					gro << species << "," << sites[i] << "," << gap_map[sites[i]] << "," << gap_map[sites[i]] << "," << file.pdb.chains[0].name << "," << aver_cor << ",0," << file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[sites[i]]].mean_pos.x << "," << file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[sites[i]]].mean_pos.y << "," << file.pdb.chains[0].aminos[file.pdb.chains[0].res_map[sites[i]]].mean_pos.z << "," << group_num << "," << group_size << endl; 
				}
			}
		}
	}
	if(site!=0){
		for(int i=0;i<sites.size();++i){
			double temp = get_distance(gap_map[site], gap_map[sites[i]], file);

			if(temp!=999.9){
				mean += temp;
				++k;
			}
		}
	}
	if(k!=0 && sites.size()>=1){
		mean /=k;
	}else{
		mean =999.9;
	}

			string temse;
			temse += species;
	for(int i=0;i<sites.size();++i){
			if(check_existance(gap_map[sites[i]], gap_map[sites[i]], file)){
				
					gro << temse << "," << sites[i] << ",-1,-1," << mean << "," << group_num+1 << endl; 
			}
	}
	return mean;
}		/* -----  end of function group_dist  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  group_dist_inter
 *  Description:  
 * =====================================================================================
 */
double group_dist_inter (vector<int>& sites, int site, Fasta_vector& file , map<int, int>& gap_map, string& species, ofstream& gro, int& group_num, unsigned int group_size)
{
	double mean=0.0;
	int k=0;



	for(int i=0;i<sites.size();++i){
		for(int j=i+1;j<sites.size();++j){	
			if(check_existance(gap_map[sites[i]], gap_map[sites[j]], file)){
				double temp = get_distance(gap_map[sites[i]], gap_map[sites[j]], file);
				
				if(temp!=999.9){
					mean += temp;
					++k;
				}
			}
		}
	}

	if(site!=0){
		for(int i=0;i<sites.size();++i){
			double temp = get_distance(gap_map[site], gap_map[sites[i]], file);

			if(temp!=999.9){
				mean += temp;
				++k;
			}
		}
	}
	if(k!=0 && sites.size()>=1){
		mean /=k;
	}else{
		mean =999.9;
	}

	for(int i=0;i<sites.size();++i){
			if(check_existance(gap_map[sites[i]], gap_map[sites[i]], file)){
				
			string temse;
			temse += species;
					gro << temse << "," << sites[i] << ",-1,-1," << mean << "," << group_num << endl; 
			}
	}
	return mean;
}		/* -----  end of function group_dist_inter  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  check_existance
 *  Description:  
 * =====================================================================================
 */
bool check_existance (int site, int site2, Fasta_vector& file)
{


	if(file.pdb.chains[0].res_map.find(site)!=file.pdb.chains[0].res_map.end() && file.pdb.chains[0].res_map.find(site2)!=file.pdb.chains[0].res_map.end()){

		return true;
	}else{
		return false;
	}


}		/* -----  end of function check_existance  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  comparison_fabs
 *  Description:  Absolute value comparison
 * =====================================================================================
 */

bool comparison_fabs(double i,double j){ 
	return (fabs(i)<fabs(j));
}


/* -----  end of function comparison_fabs  ----- */
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  SD_vf
 *  Description:  Standard deviation of a vector of doubles
 * =====================================================================================
 */


double SD_vf(vector<double>& array, double mean){
	double total=0.0;
	int size = array.size();

	for(int i=0;i<size;++i){

		total += pow(((double)array[i] - mean),2);

	}
	total =total/size;

	return sqrt(total);
}

/* -----  end of function SD_vf  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Diff_aa_column
 *  Description:  find differences in the columns
 * =====================================================================================
 */


int Diff_aa_column(vector<string>& sequences, vector<int>& array){;

	int seq_number = sequences.size();
	int length_sequences = sequences[0].size();
	array.assign(length_sequences, 0);

	for(int i=0;i<length_sequences;++i){
		map<char, int> diffmap;
		for(int j=0;j < seq_number;++j){
			/*	if(i==0 && statl==1){
				cerr << "daa=" << sequences[j][i] << endl;
				}*/
			map<char, int>::iterator it;
			it = diffmap.find((char)sequences[j][i]);
			if(it==diffmap.end()){
			diffmap[(char)sequences[j][i]]=1;
			}else{
				++diffmap[(char)sequences[j][i]];
			}
			

			if(diffmap.size()>1){
			int total=0;
			map<char, int>::iterator mit;
				for(mit=diffmap.begin();mit!=diffmap.end();++mit){
					if(mit->second>=2){
						++total;
					}
				}
				if(total>=2){
				array[i]=1;
				break;
				}
	
			}
		}//end while j
	}//end while
	return 0;
}


/* -----  end of function Diff_aa_column  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  D_average
 *  Description:  Gets the average of each row of D
 * =====================================================================================
 */


int D_average(vector< vector<double> >& D, vector<double>& D_aver){
	int i=0, j=0, totals;

	totals = D[0].size();
	int size_D = D.size();

	D_aver.assign(size_D, 0.0);

	while(i<size_D-1){
		j=0;
		while(j<=totals){
			D_aver[i] += D[i][j];
			++j;

		}
		D_aver[i]/=totals;
		++i;
	}
	return 0;
}

/* -----  end of function D_average  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  average_vi
 *  Description:  Gives average of a vector of ints, returns double
 * =====================================================================================
 */

double average_vi(vector<int>& array){
	double total;
	int size = array.size();

	total=0.0;
	for(int i=0; i<size;++i){

		total += array[i];
	}
	total /= size;
	return total; 
}
/* -----  end of function average_vi  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getAncestralSequences
 *  Description:  
 * =====================================================================================
 */



int getAncestralSequences(int args, char ** argv, TreeTemplate<Node>* mytree, vector<string>& sequences, vector<string>& names, map<int, string>& mapinner)
{

	try {

		BppApplication bppancestor(args, argv, "BppAncestor");
		bppancestor.startTimer();

		Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppancestor.getParams(), "", false);

		//		VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppancestor.getParams());
		VectorSiteContainer *allSites = new VectorSiteContainer(alphabet);


		SequenceContainer *myseq = SequenceContainerTools::createContainerOfSpecifiedSize(alphabet, sequences.size());

		//SequenceContainer *myseq = new SequenceContainer;

		for(int i=0;i<sequences.size();++i){
			Sequence temseq(names[i], sequences[i], alphabet);
			myseq->addSequence(temseq, true);
			allSites->addSequence(temseq, alphabet);
		}	
		delete myseq;
		//cerr << "FIRST" << endl;


		VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, bppancestor.getParams());
		delete allSites;
		//		ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
		//		ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));

		// Get the initial tree
		//Tree* tree = PhylogeneticsApplicationTools::getTree(bppancestor.getParams());
		Tree* tree = mytree->clone();

		//		ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));

		string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppancestor.getParams(), false, false);
		if(treeWIdPath != "none")
		{
			TreeTemplate<Node> ttree(*tree);
			vector<Node *> nodes = ttree.getNodes();
			for(unsigned int i = 0; i < nodes.size(); i++)
			{
				if(nodes[i]->isLeaf())
					nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
				else
					nodes[i]->setBranchProperty("NodeId", BppString(TextTools::toString(nodes[i]->getId())));
			}
			Newick treeWriter;
			treeWriter.enableExtendedBootstrapProperty("NodeId");
			ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
			treeWriter.write(ttree, treeWIdPath);
			delete tree;
			cout << "BppML's done." << endl;
			exit(0);
		}

		DRTreeLikelihood *tl;
		string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppancestor.getParams(), "no", "", true, false);
		ApplicationTools::displayResult("Heterogeneous model", nhOpt);

		SubstitutionModel    *model    = 0;
		SubstitutionModelSet *modelSet = 0;
		DiscreteDistribution *rDist    = 0;
		unsigned int nbStates;

		if(nhOpt == "no")
		{  
			model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppancestor.getParams());
			if(model->getNumberOfStates() > model->getAlphabet()->getSize())
			{
				//Markov-modulated Markov model!
				rDist = new ConstantDistribution(1.);
			}
			else
			{
				rDist = PhylogeneticsApplicationTools::getRateDistribution(bppancestor.getParams());
			}
			tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true);
			nbStates = model->getNumberOfStates();
		}
		else if(nhOpt == "one_per_branch")
		{
			model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppancestor.getParams());
			if(model->getNumberOfStates() > model->getAlphabet()->getSize())
			{
				//Markov-modulated Markov model!
				rDist = new ConstantDistribution(1.);
			}
			else
			{
				rDist = PhylogeneticsApplicationTools::getRateDistribution(bppancestor.getParams());
			}
			vector<double> rateFreqs;
			if(model->getNumberOfStates() != alphabet->getSize())
			{
				//Markov-Modulated Markov Model...
				unsigned int n =(unsigned int)(model->getNumberOfStates() / alphabet->getSize());
				rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
				// we should assume a rate distribution for the root also!!!  
			}
			FrequenciesSet * rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, sites, bppancestor.getParams(), rateFreqs);
			vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", bppancestor.getParams(), ',', "");
			modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters); 
			model = 0;
			tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
			nbStates = modelSet->getNumberOfStates();
		}
		else if(nhOpt == "general")
		{
			modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, sites, bppancestor.getParams());
			if(modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
			{
				//Markov-modulated Markov model!
				rDist = new ConstantDistribution(1.);
			}
			else
			{
				rDist = PhylogeneticsApplicationTools::getRateDistribution(bppancestor.getParams());
			}
			tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
			nbStates = modelSet->getNumberOfStates();
		}
		else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);

		tl->initialize();

		delete tree;

		double logL = tl->getValue();
		if(isinf(logL))
		{
			// This may be due to null branch lengths, leading to null likelihood!
			ApplicationTools::displayWarning("!!! Warning!!! Likelihood is zero.");
			ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
			ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
			ParameterList pl = tl->getBranchLengthsParameters();
			for(unsigned int i = 0; i < pl.size(); i++)
			{
				if(pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
			}
			tl->matchParametersValues(pl);
			logL = tl->getValue();
		}
		if(isinf(logL))
		{
			ApplicationTools::displayError("!!! Unexpected likelihood == 0.");
			ApplicationTools::displayError("!!! Looking at each site:");
			for(unsigned int i = 0; i < sites->getNumberOfSites(); i++)
			{
				(*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
			}
			ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
			exit(-1);
		}
		tree = new TreeTemplate<Node>(tl->getTree());

		// Write parameters to screen:
		ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
		ParameterList parameters = tl->getSubstitutionModelParameters();
		for(unsigned int i = 0; i < parameters.size(); i++)
		{
			ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
		}
		parameters = tl->getRateDistributionParameters();
		for(unsigned int i = 0; i < parameters.size(); i++)
		{
			ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
		}

		// Getting posterior rate class distribution:
		DiscreteDistribution* prDist = RASTools::getPosteriorRateDistribution(*tl);
		ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
		if (ApplicationTools::message) prDist->print(*ApplicationTools::message);
		ApplicationTools::displayMessage("\n");
		delete prDist;

		// Reconstruct ancestral sequences:
		string reconstruction = ApplicationTools::getStringParameter("asr.method", bppancestor.getParams(), "marginal", "", true, false);
		ApplicationTools::displayResult("Ancestral state reconstruction method", reconstruction);
		bool probs = false;

		AncestralStateReconstruction *asr = NULL;
		bool probMethod = false;
		if(reconstruction == "marginal")
		{
			asr = new MarginalAncestralStateReconstruction(tl);
			probMethod = true;
		}
		else
			throw Exception("Unknown ancestral state reconstruction method: " + reconstruction);

		if(probMethod)
		{
			probs = ApplicationTools::getBooleanParameter("asr.probabilities", bppancestor.getParams(), false, "", true, false);
			ApplicationTools::displayResult("Output probabilities", probs ? "yes" : "no");
		}

		// Write infos to file:
		string outputFile = ApplicationTools::getAFilePath("output.sites.file", bppancestor.getParams(), false, false);
		if(outputFile != "none")
		{
			ApplicationTools::displayResult("Output file for sites", outputFile);
			ofstream out(outputFile.c_str(), ios::out);
			TreeTemplate<Node> ttree(*tree);
			vector<Node *> nodes = ttree.getInnerNodes();
			unsigned int nbNodes = nodes.size();

			// Get the rate class with maximum posterior probability:
			vector<unsigned int> classes = tl->getRateClassWithMaxPostProbOfEachSite();
			// Get the posterior rate, i.e. rate averaged over all posterior probabilities:
			Vdouble rates = tl->getPosteriorRateOfEachSite();
			// Get the ancestral sequences:
			vector<Sequence *> sequences(nbNodes);
			vector<VVdouble *> probabilities(nbNodes);

			vector<string> colNames;
			colNames.push_back("Sites");
			colNames.push_back("is.complete");
			colNames.push_back("is.constant");
			colNames.push_back("lnL");
			colNames.push_back("rc");
			colNames.push_back("pr");
			for(unsigned int i = 0; i < nbNodes; i++)
			{
				Node *node = nodes[i];
				colNames.push_back("max." + TextTools::toString(node->getId()));
				if(probs)
				{
					probabilities[i] = new VVdouble();
					//The cast will have to be updated when more probabilistic method will be available:
					sequences[i] = dynamic_cast<MarginalAncestralStateReconstruction *>(asr)->getAncestralSequenceForNode(node->getId(), probabilities[i], false);

					for(unsigned int j = 0; j < nbStates; j++)
					{
						colNames.push_back("prob." + TextTools::toString(node->getId()) + "." + alphabet->intToChar((int)j));
					}
				}
				else
				{
					sequences[i] = asr->getAncestralSequenceForNode(node->getId());
					//cerr << "sequences[" << i << "]=" << sequences[i] << endl;
				}
			}

			//Now fill the table:
			vector<string> row(colNames.size());
			DataTable* infos = new DataTable(colNames);

			for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
			{
				double lnL = tl->getLogLikelihoodForASite(i);
				const Site* currentSite = &sites->getSite(i);
				int currentSitePosition = currentSite->getPosition();
				int isCompl = (SiteTools::isComplete(* currentSite) ? 1 : 0);
				int isConst = (SiteTools::isConstant(* currentSite) ? 1 : 0);
				row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
				row[1] = TextTools::toString(isCompl);
				row[2] = TextTools::toString(isConst);
				row[3] = TextTools::toString(lnL);
				row[4] = TextTools::toString(classes[i]);
				row[5] = TextTools::toString(rates[i]);

				unsigned int k = 6;
				for(unsigned int j = 0; j < nbNodes; j++)
				{
					row[k] = sequences[j]->getChar(i);
					k++;
					if(probs)
					{
						for(unsigned int l = 0; l < nbStates; l++)
						{
							row[k] = TextTools::toString((*probabilities[j])[i][l]);
							k++;
						}
					}
				}

				infos->addRow(row);
			}

			DataTable::write(*infos, out, "\t");

			delete infos;
		}


		outputFile = ApplicationTools::getAFilePath("output.nodes.file", bppancestor.getParams(), false, false);
		/*if (outputFile != "none")
		  {
		  ApplicationTools::displayResult("Output file for nodes", outputFile);
		  ofstream out(outputFile.c_str(), ios::out);

		  map<int, vector<double> > frequencies;
		  TreeLikelihoodTools::getAncestralFrequencies(*tl, frequencies, false);

		  vector<string> colNames;
		  colNames.push_back("Nodes");
		  for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
		  colNames.push_back("exp" + TextTools::toString(i));
		  for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
		  colNames.push_back("eb" + TextTools::toString(i));

		//Now fill the table:
		vector<string> row(colNames.size());
		DataTable* infos = new DataTable(colNames);

		for (map<int, vector<double> >::iterator it = frequencies.begin(); it != frequencies.end(); it++)
		{
		row[0] = TextTools::toString(it->first);
		Vdouble ebFreqs = DRTreeLikelihoodTools::getPosteriorStateFrequencies(*tl, it->first);
		for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
		{
		row[i + 1] = TextTools::toString(it->second[i]);
		}
		for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
		{
		row[i + tl->getNumberOfStates() + 1] = TextTools::toString(ebFreqs[i]);
		}
		infos->addRow(row);
		}

		DataTable::write(*infos, out, "\t");

		delete infos;
		}
		*/


		SequenceContainer *asSites = 0;
		if(probMethod)
		{
			bool sample = ApplicationTools::getBooleanParameter("asr.sample", bppancestor.getParams(), false, "", true, false);
			ApplicationTools::displayResult("Sample from posterior distribution", sample ? "yes" : "no");
			if(sample)
			{
				unsigned int nbSamples = ApplicationTools::getParameter<unsigned int>("asr.sample.number", bppancestor.getParams(), 1, "", true, false);
				asSites = new AlignedSequenceContainer(alphabet);
				for(unsigned int i = 0; i < nbSamples; i++)
				{
					ApplicationTools::displayGauge(i, nbSamples-1, '=');
					SequenceContainer *sampleSites = dynamic_cast<MarginalAncestralStateReconstruction *>(asr)->getAncestralSequences(true);
					vector<string> names = sampleSites->getSequencesNames();
					for(unsigned int j = 0; j < names.size(); j++)
						names[j] += "_" + TextTools::toString(i+1);
					sampleSites->setSequencesNames(names, false);
					SequenceContainerTools::append(*asSites, *sampleSites);
					delete sampleSites;
				}
				ApplicationTools::message->endLine();
			}
			else
			{
				asSites = asr->getAncestralSequences();
			}
		}
		else
		{
			asSites = asr->getAncestralSequences();
		}

		vector<int> inner;
		TreeTemplate<Node> trel(*tree);
		inner = TreeTemplateTools::getInnerNodesId(*trel.getNode(trel.getRootId()));

		for(int i=0;i<inner.size();++i){
			Sequence  *my = asr->getAncestralSequenceForNode(inner[i]);
			//			cerr << "inner[" << i << "]=" << inner[i] << my->toString()  << endl;
			mapinner[inner[i]] = my->toString(); 
			delete my;
		}
		//	cerr << endl;

		//		SequenceApplicationTools::writeSequenceFile(*asSites, bppancestor.getParams());

		delete asSites;

		delete asr;
		delete alphabet;
		delete sites;
		if(model)    delete model;
		if(modelSet) delete modelSet;
		delete rDist;
		delete tl;
		delete tree;
		bppancestor.done();

	}
	catch(exception & e)
	{
		cout << e.what() << endl;
		return 1;
	}

	return 0;
}

/* -----  end of function getAncestralSequences  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Estimate_D
 *  Description:  Estimation of D
 * =====================================================================================
 */


void Estimate_D(vector<string>& seqes, vector<string>& names, vector< vector<double> >& D, vector< vector<double> >& D_corr, vector<double>& rel_dist, int correction, Blosum& Blos62, vector<int>& aa_differences, TreeTemplate<Node>* tree){

	int number_thetas = 0, total_theta;
	double Theta_aver;



	char **arg;


	int ar = 8;

	arg = (char**)malloc(ar*sizeof(char *));

	for(int i=0;i<ar;++i){
		arg[i] = (char*)malloc(256*sizeof(char));
	}

	sprintf(arg[0], "%s", "bppancestor");
	sprintf(arg[1], "%s", "input.sequence.file=YAL017W.seq.alnc.aln");
	sprintf(arg[2], "%s", "alphabet=Protein");
	sprintf(arg[3], "%s", "input.sequence.sites_to_use=all");
	sprintf(arg[4], "%s", "input.sequence.max_gap_allowed=100%");
	sprintf(arg[5], "%s", "input.tree.file=yeast_tree.tre");
	sprintf(arg[6], "%s", "model=JTT92");
	sprintf(arg[7], "%s", "output.sequence.file=test2");

	map<int, string> mapinner;
	getAncestralSequences(ar, arg, tree, seqes, names, mapinner);

	for(int i=0;i<ar;++i){
		free(arg[i]);
	}
	free(arg);



	int num_seqs = seqes.size();
	int length_seqs = seqes[0].size();
	total_theta = num_seqs*(num_seqs-1);
	//	int number_comp = (num_seqs*(num_seqs -1))/2;

	vector<int> nodes = tree->getNodesId();
	int number_comp  = nodes.size()-1;

	vector<double> tempD(number_comp, 0.0);
	vector<int> leaves = tree->getLeavesId();
	map<string, int> names_id;

	for(int i=0;i<leaves.size();++i){
		names_id[tree->getNodeName(leaves[i])] = leaves[i];
	}

	for(int i=0;i<length_seqs;++i){
		D.push_back(tempD);

		int k=0;
		vector<int> theta;
		vector<int> fathers;
		int test=1;
		for(int j=0;j<names.size();++j){

			int nod = names_id[names[j]]; 
			int darth = tree->getFatherId(nod);
			if(darth!=tree->getRootId()){
				//		fathers.push_back(darth);
				string temp, temp2;
				temp += mapinner[darth][i];
				temp2 += seqes[j][i];
				//			int test = Blos62.values[Blos62.index[temp2]][Blos62.index[temp]];
				theta.push_back(Blos62.values[Blos62.index[temp2]][Blos62.index[temp]]);
				//				theta.push_back(test);
			}
		}

		map<int, string>::iterator mit;

		for(mit=mapinner.begin();mit!=mapinner.end();++mit){

			if(mit->first!=tree->getRootId()){
				string temp, temp2;
				int da = tree->getFatherId(mit->first);	
				temp += mapinner[da][i];
				temp2  += mapinner[mit->first][i];
				theta.push_back(Blos62.values[Blos62.index[temp2]][Blos62.index[temp]]);	
			}

		}


		double Theta_av = average_vi(theta);

		for(int j=0;j<theta.size();++j){
			D[i][j] = pow((theta[j] - Theta_av),2);
			//				D[i][j] = theta[j];
		}

	}

}

/* -----  end of function Estimate_D  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Correlation
 *  Description:  Pearsons R
 * =====================================================================================
 */

double Correlation(vector<double>& sample_a, vector<double>& sample_b){
	double mean_a=0.0, mean_b=0.0, nominator=0.0, X=0.0, Y=0.0, correl=0.0, test;


	mean_a = average_vec<double>(sample_a);
	mean_b = average_vec<double>(sample_b);

	for(int i=0;i<sample_a.size();++i){
		nominator += (sample_a[i] - mean_a) * (sample_b[i] - mean_b);
		//		cerr << "sample_b[" << i << "]=" << sample_b[i] << endl;
		Y += (pow((sample_b[i] - mean_b), 2));
		//	}
		/*else{
		  nominator+=0;
		  Y+=0;
		  }*/
		X += (pow((sample_a[i] - mean_a),2));
	}

	//cerr << "mean_B=" << mean_b << endl;
	//cerr << "X=" << X << "\tY=" << Y << endl;
	/*/end while i*/
	test = sqrt(X*Y);
	if(test!=0.0){
		correl = (nominator)/test;
	}/*else
	   {
	   correl = 0.0;
	   }*/
	//cerr << "correl=" << correl << endl;
	return correl;
}

/* -----  end of function Correlation  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  intra
 *  Description:  Intra-molecular co-evolution
 * =====================================================================================
 */


double intra(vector< vector<double> >& D, vector< vector<double> >& D_correct, vector<double>& Correl, int information13, int length_seq, int simulate){
	int j=1, m=0, n=0, number_correlations=0, num =0, size_D_column, size_D_column2, total_number_comp;
	double percent;

	size_D_column = D[0].size();

	total_number_comp =(( (length_seq) - 1) * length_seq)/2;
	vector<double> D_column(size_D_column, 0.0), D_column2(size_D_column, 0.0);



	/* change these loops so that we have the whole i-loop twice within each of if and else, profiling!! */
	for(int i=0;i<D.size()-1;++i){
		m = 0;

		if(information13 == 0){
			while(m <= size_D_column - 1){
				D_column[m] = D[i][m];
				m++;
			}

		}else{

			while(m <= size_D_column - 1){
				D_column[m] = D_correct[i][m];
				m++;
			}
		}

		for(j=i+1;j< D.size();++j){
			n = 0;

			if(information13 == 0){
				while(n <= size_D_column - 1){
					D_column2[n] = D[j][n];
					n++;
				}
			}else{
				while(n <= size_D_column - 1){
					D_column2[n] = D_correct[j][n];
					n++;
				}
			}

			if(simulate==1){
				Correl.push_back(fabs(Correlation(D_column, D_column2)));
			}else{
				Correl.push_back(Correlation(D_column, D_column2));
			}

			++number_correlations;
		}
		percent = ((double)i)/ ((double)(length_seq-1)) * 100;
		j = i + 1;
	}

	return number_correlations;
}


/* -----  end of function intra  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  create_input_tree
 *  Description:  makes an input tree which can be scored
 * =====================================================================================
 */

TreeTemplate<Node> *create_input_tree(vector<string>& seq_names, vector< string >& sequences){
	FILE *stream;
	stream = freopen("/dev/null", "w", stdout);

	std::auto_ptr<DistanceMatrix> DS;
	vector<string> names;

	int tot=0;
	names.resize(seq_names.size());
	for(int i=0;i<(signed)seq_names.size();++i){
		names[i] = seq_names[i];
		++tot;
	}

	DS=ScoreDist(names, sequences);
	/*if(DS==NULL){
		return NULL;
	}*/

	AgglomerativeDistanceMethod * distMethod = NULL;
	BioNJ * bionj = new BioNJ();
	bionj->outputPositiveLengths(true);
	distMethod = bionj;

	bionj->setDistanceMatrix(*DS);
	bionj->computeTree(true);
	//if(DS!=NULL)
	//	delete DS;
	TreeTemplate<Node> * tree2 = dynamic_cast<TreeTemplate<Node> *>(bionj->getTree());

	delete bionj;
	return tree2;
}


void read_hydros(char *filename){
	fstream hyd;

	hyd.open(filename);





}

int print_splash(string filename){


	if(filename=="stdout"){

		fprintf(stdout, "\n\t\t*******************************************************************************************\n\t\t");
		fprintf(stdout, "* CAPS: Co-Evolution Analysis using Protein Sequences                                     *\n\t\t");
		fprintf(stdout, "* Author: Brian E. Caffrey								  *\n\t\t");
		fprintf(stdout, "* Evolutionary Genetics and Bioinformatics Laboratory					  *\n\t\t");
		fprintf(stdout, "* Smurfit Institute of Genetics								  *\n\t\t");
		fprintf(stdout, "* University of Dublin, Trinity College							  *\n\t\t");
		fprintf(stdout, "* Mathematical Model: Caffrey, Williams, Hokamp and Fares, 2011		  			  *\n\t\t");
		fprintf(stdout, "*******************************************************************************************\n");
	}else{

		ofstream OUTPUT(filename.c_str());
		OUTPUT << "\n\t\t*******************************************************************************************\n\t\t";
		OUTPUT << "* CAPS: Co-Evolution Analysis using Protein Sequences                                     *\n\t\t";
		OUTPUT << "* Author: Brian E. Caffrey									  *\n\t\t";
		OUTPUT << "* Evolutionary Genetics and Bioinformatics Laboratory					  *\n\t\t";
		OUTPUT << "* Smurfit Institute of Genetics								  *\n\t\t";
		OUTPUT << "* University of Dublin, Trinity College							  *\n\t\t";
		OUTPUT << "* Mathematical Model: Caffrey, Williams, Hokamp and Fares, 2011  					*\n\t\t";
		OUTPUT << "*******************************************************************************************\n";
	}

	return 0;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Remove_names
 *  Description:  
 * =====================================================================================
 */
TreeTemplate<Node>* Remove_names (vector<string>& names, TreeTemplate<Node>* tree){

	vector<int> nodes = tree->getNodesId();

	for(unsigned int i=0;i<nodes.size();++i){
		tree->deleteDistanceToFather(nodes[i]);
	}


	vector<string> treenames = tree->getLeavesNames();
	vector<int> indices;

	/*get indices we don't want*/

	for(unsigned int i=0;i<treenames.size();++i){
		int found =0;
		for(unsigned int j=0;j<names.size();++j){

			if(names[j]==treenames[i]){
				found =1;
				break;
			}
		}

		if(found==0){
			indices.push_back(i);
		}

	}


	string treestring = TreeTemplateTools::treeToParenthesis(*tree, false);

	for(unsigned int i=0;i<indices.size();++i){

		size_t found;
		found = treestring.find(treenames[indices[i]]);

		if(found!=string::npos){
			int start = int(found);
			int end = int(found) + treenames[indices[i]].size();

			string left = whats_left(start, treestring);
			string right = whats_right(end, treestring);

			if(left=="comma" && right=="bracket"){
				int pos;
				find_left_bracket(pos, start, treestring);
				treestring.erase(start+1, end-start);
				treestring.erase(pos, 1);

			}else if(left=="bracket" && right=="comma"){
				int pos;
				find_right_bracket(pos, end, treestring);
				treestring.erase(pos, 1);
				treestring.erase(start, end-start+1);

			}else if(left=="bracket" && right=="bracket"){
				int pos;
				find_comma(pos, start, end, treestring);
				if(pos<start){
					treestring.erase(start, end-start);
					treestring.erase(pos, 1);
				}else{
					treestring.erase(pos, 1);
					treestring.erase(start, end-start);
				}
			}

		}

	}

	TreeTemplate<Node> *newtree;
	stringstream strstr(treestring);

	Newick * newickReader = new Newick(false);
	newtree = newickReader->read(strstr);
	//free(strstr);
	delete newickReader;

	return newtree;
}		/* -----  end of function Remove_names  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_comma
 *  Description:  
 * =====================================================================================
 */
void find_comma (int& pos, int start, int end, string tree){

	int left=0, right=0;

	for(int i=start;i>=0;--i){
		if(tree[i]==','){
			left = start - i;
			break;
		}
	}

	for(unsigned int i=end;i<tree.size();++i){
		if(tree[i]==','){
			right =  i - end;
			break;
		}
	}

	if(left<right){
		pos = start - left;
	}else{
		pos = end + right;
	}

}		/* -----  end of function find_comma  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_right_bracket
 *  Description:  
 * =====================================================================================
 */
void find_right_bracket (int& pos, int end, string tree){

	int right=0, left=0;

	for(unsigned int i=end;i<tree.size();++i){

		if(tree[i]==')'){
			++right;
		}
		if(tree[i]=='('){
			++left;
		}
		if(right>left){
			pos = i;
			return;
		}
	}


}		/* -----  end of function find_left_bracket  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_left_bracket
 *  Description:  
 * =====================================================================================
 */
void find_left_bracket (int& pos, int start, string tree){

	//cerr << "HI" << endl;
	//exit(-1);

	int right=0, left=0;

	for(int i=start;i>=0;--i){

		if(tree[i]==')'){
			++right;
		}
		if(tree[i]=='('){
			++left;
		}
		if(left>right){
			pos = i;
			return;
		}
	}


}		/* -----  end of function find_left_bracket  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  whats_right
 *  Description:  
 * =====================================================================================
 */
string whats_right (int& end, string treestring){

	for(unsigned int i=end;i<treestring.size();++i){

		if(treestring[i]==','){
			end = i;
			return (string)"comma";
		}else if(treestring[i]==')'){
			end = i;
			return (string)"bracket";
		}

	}

	return "error";
}		/* -----  end of function whats_right  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  whats_left
 *  Description:  what is left of the name in the Newick string
 * =====================================================================================
 */
string whats_left (int& start, string tree){

	for(int i=start;i>=0;--i){
		if(tree[i]==','){
			start = i;
			return (string)"comma";
		}else if(tree[i]=='('){
			start = i;
			return (string)"bracket";
		}
	}
	return "error";
}		/* -----  end of function whats_left  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Put_distances_on_tree
 *  Description:  
 * =====================================================================================
 */


vector<Node *> Put_distances_on_tree(TreeTemplate<Node> *tree, std::auto_ptr<DistanceMatrix> DS ){	

	vector<int> nodes = tree->getNodesId();
	tree->setBranchLengths(0.0);

	//vector<string> nams;
	//	nams = DS->getNames();
	map<string, int> present;
	vector<int> tips;
	vector<Node *> pres;
	vector<int> fathers;
	vector<Node *> fath;

	vector<string> names = tree->getLeavesNames();


	/* make nested nodes using the vector of ints of leaves */
	vector<int> isLeaf;
	isLeaf = tree->getLeavesId();





	/*put the distances on the tips */
	vector<int> temp_distances;
	for(unsigned int j=0;j<isLeaf.size()-2;++j){
		pres.push_back(tree->getNode(isLeaf[j]));
		double dist1, dist2, dist3, dist4;

		dist1 = (*DS)(tree->getNodeName(isLeaf[j]), tree->getNodeName(isLeaf[j+1]));
		dist2 = (*DS)(tree->getNodeName(isLeaf[j]), tree->getNodeName(isLeaf[j+2]));
		dist3 = (*DS)(tree->getNodeName(isLeaf[j+1]), tree->getNodeName(isLeaf[j+2]));
		dist4 = (dist1+dist2-dist3)/2;
		tree->setDistanceToFather(isLeaf[j], dist4);
		temp_distances.push_back(dist4);
		dist4 = (dist1-dist2+dist3)/2;
		tree->setDistanceToFather(isLeaf[j+1], dist4);
		temp_distances.push_back(dist4);
		dist4 = (-dist1+dist2+dist3)/2;
		tree->setDistanceToFather(isLeaf[j+2], dist4);
		temp_distances.push_back(dist4);
	}
	pres.push_back(tree->getNode(isLeaf[isLeaf.size()-2]));
	pres.push_back(tree->getNode(isLeaf[isLeaf.size()-1]));

	vector<int> inner;
	vector<Node *> innod;

	innod = tree->getInnerNodes();
	inner = TreeTemplateTools::getInnerNodesId(*tree->getRootNode());

	map<int, int> idmap;
	map<Node *, int> nodmap;
	int j=0;
	for(unsigned int i=0;i<inner.size();++i){
		if(inner[i]!=tree->getRootId()){

			int da = tree->getFatherId(inner[i]);
			vector<int> sol = tree->getSonsId(da);
			if(sol.size()>1){
				idmap[inner[i]] = j;
				++j;
			}else{
				tree->setDistanceToFather(sol[0], 0.0);
			}
		}
	}

	vector<int> sonny = tree->getSonsId(tree->getRootId());
	vector< vector<double> > pol;
	vector< double> temper(idmap.size(), 0.0);
	vector< double> ans(idmap.size(), 0.0);
	for(unsigned int i=0;i<idmap.size();++i){
		pol.push_back(temper);
	}


	map<int, int> combo;
	int k=0;

	map<int, int>::iterator mi;

	for(mi=idmap.begin();mi!=idmap.end();++mi){

		get_unique_pair(combo, mi->first, tree, idmap, pol, k, ans, DS);
		++k;
	}

	gsl_matrix *A; 
	gsl_vector *x, *b;

	A= gsl_matrix_alloc(pol.size(), pol[0].size());
	b= gsl_vector_alloc(pol[0].size());
	x= gsl_vector_alloc(pol[0].size());


	for(unsigned int i=0;i<pol.size();++i){
		for(unsigned int j=0;j<pol[i].size();++j){
			gsl_matrix_set(A, i, j, pol[i][j]);
		}
		gsl_vector_set(b, i, ans[i]);
	}

	gsl_linalg_HH_solve (A, b, x);

	for(unsigned int i=0;i<pol.size();++i){

		if(gsl_vector_get(x, i)<=0){
			gsl_vector_set(x, i, 0.01);
		}


		if((double)gsl_vector_get(x,i)<=0.0){
			tree->setDistanceToFather(inner[i], gsl_vector_get(x, i));	
		}else{
			tree->setDistanceToFather(inner[i], 0.001);	
		}

	}


	string tempp = TreeTemplateTools::treeToParenthesis(*tree, false);




	for(unsigned int i=0;i<nodes.size();++i){

		if(tree->hasDistanceToFather(nodes[i])==false){
			if(nodes[i]==tree->getRootId()){
				tree->setDistanceToFather(nodes[i], 0.0);
			}else{
				cerr << "Error setting distances on tree!" << endl;
				exit(-1);
			}
		}

	}

	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_matrix_free(A);


	return pres;
}
/* -----  end of function Put_distances_on_tree  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  get_unique_pair
 *  Description:  
 * =====================================================================================
 */
void get_unique_pair (map<int, int>& combo, int node, TreeTemplate<Node>* tree, map<int, int>& idmap, vector< vector<double> >& pol, int l, vector<double>& ans, std::auto_ptr<DistanceMatrix> DS ){

	Node *subtreeRoot = TreeTemplateTools::cloneSubtree<Node>(*tree->getNode(node));
	TreeTemplate<Node> subtree(subtreeRoot); //Create a new (independent) tree from the subtree
	vector<int> leavesid = subtree.getLeavesId();


	map<int, int>::iterator myit;
	myit = combo.find(leavesid[0]);

	if(myit!=combo.end() && myit->second!=leavesid[leavesid.size()-1]){
		combo[myit->first]=leavesid[leavesid.size()-1];
	}

	int g=0;

	if(tree->hasNodeName(leavesid[0])==false && leavesid.size()==1){
		pol.erase(pol.begin()+l, pol.begin()+l+1);
		return;
	}

	while(tree->hasNodeName(leavesid[g])!=true){
		++g;
	}


	int h = leavesid.size()-1;

	while(tree->hasNodeName(leavesid[h])!=true){
		--h;
		if(h==0){
			return;
		}
	}

	int dar1, dar2;

	dar1=tree->getFatherId(leavesid[g]);
	map<int, int>::iterator mit;
	mit = idmap.find(dar1);
	if(dar1!=node)
		pol[l][mit->second]=1.0;
	while(dar1!=node){
		dar1=tree->getFatherId(dar1);
		if(dar1!=node)
			pol[l][mit->second]=1.0;
	}



	dar2=tree->getFatherId(leavesid[h]);
	mit = idmap.find(dar2);
	if(dar2!=node)
		pol[l][mit->second]=1.0;
	while(dar2!=node){
		dar2=tree->getFatherId(dar2);
		if(dar1!=node)
			pol[l][mit->second]=1.0;
	}


	pol[l][idmap[node]]=1.0;
	ans[l] -= tree->getDistanceToFather(leavesid[g]);
	ans[l] -= tree->getDistanceToFather(leavesid[h]);


	ans[l] += (*DS)(tree->getNodeName(leavesid[g]), tree->getNodeName(leavesid[h]));
}		/* -----  end of function get_unique_pair  ----- */


