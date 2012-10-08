#include"Blosum.h"
#include<algorithm>

typedef struct{

int pos[2],real[2];

}site;


typedef struct{

	int size;
	char **nuc, **amino;

}NUC_AA;


typedef struct{
	int size, pos, real;
	char string[100000];
}interstring;

template<class B>
double average_vec(vector<B>& array){
	double total=0.0;

	for(unsigned int i=0;i<array.size();++i){

		total += array[i];
	}
	total /= (double)array.size();
	return total; 
}


class Nod{
	public:
		map<int, vector<int> > node_connect;
};



class Network{

	public:
		void add_connection(int& i, int& j);
		vector<Nod> nodes;	
		map<int, int> node_map;
		void print_network(ofstream& file);	


};




class hydro{
	public:
		void set_hyd ();
		map<char, int> values;
};

 
class MW{
	public:
		void set_mw ();
		map<char, int> values;
};


string assess_hyd(vector<string>& sequences, int col1, int col2);
void Hydro (vector< vector<int> >& pairs, vector<string>& sequences, ofstream& output);
string whats_left (int& start, string tree);
string whats_right (int& end, string treestring);
	void find_right_bracket (int& pos, int end, string tree);
	void find_left_bracket (int& pos, int start, string tree);
void find_comma (int& pos, int start, int end, string tree);
TreeTemplate<Node>* Remove_names (vector<string>& names, TreeTemplate<Node>* tree);
//vector<Node *> Put_distances_on_tree(TreeTemplate<Node> *tree, DistanceMatrix *DS );	
vector<Node *> Put_distances_on_tree(TreeTemplate<Node> *tree, std::auto_ptr<DistanceMatrix> DS );	
void get_unique_pair (map<int, int>& combo, int node, TreeTemplate<Node>* tree, map<int, int>& idmap, vector< vector<double> >& pol, int l, vector<double>& ans, std::auto_ptr<DistanceMatrix> DS );

hydro sethyd();
	bool check_existance (int site, int site2, Fasta_vector& file);
double group_dist_intra (vector<int>& sites, int site, Fasta_vector& file , map<int, int>& gap_map, string& species, ofstream& gro, int& group_num, unsigned int group_size);
double group_dist_inter (vector<int>& sites, int site, Fasta_vector& file , map<int, int>& gap_map, string& species, ofstream& gro, int& group_num, unsigned int group_size);
double  get_distance (int pos1, int pos2, Fasta_vector& file);
void group_inter (vector< vector<int> >& pairs, map<int, int >& gap_map1, map<int, int>& gap_map2, ofstream& output, vector<string>& sequences, vector<string>& sequences2, Fasta_vector& file1, Fasta_vector& file2, string& ref_spec, string& ref_spec2, string& outfile, string& files1, string& files2);
void group_intra (vector< vector<int> >& pairs, ofstream& output, Fasta_vector& file, map<int, int>& gap_map, string& ref_spec, string& filei);
int getAncestralSequences(int args, char ** argv, TreeTemplate<Node>* mytree, vector<string>& sequences, vector<string>& names, map<int, string>& mappinner);
void Estimate_D(vector<string>& seqes, vector<string>& names, vector< vector<double> >& D, vector< vector<double> >& D_corr, vector<double>& rel_dist, int correction, Blosum& Blos62, vector<int>& aa_differences, TreeTemplate<Node>* tree);
	bool check_gaps (vector<string>& seq, vector<string>& seq2, int i, int j);
int Chi_squared (int num_pairs, int num_correlations, double background, double& P_val);
int print_inter(vector<double>& Correl, double threshold, ofstream& output, vector<string>& sequences, vector<string>& sequences2, vector< vector<double> >& D, vector< vector<double> >& D2, vector<int>& diff1, vector<int>& diff2, Fasta_vector& file1, Fasta_vector& file2, string& ref_spec1, string& ref_spec2, string& outfile, string& files1, string& files2, int ref, double& thresholdR);
	void non_over_intra (vector< vector<int> >& pairs, ofstream& output, Fasta_vector& file, map<int, int>& gap_map,      string& ref_spec, string& filei);
	bool not_subset (vector< vector<int> >& subs, vector<int>& set );
void non_over_inter (vector< vector<int> >& pairs, map<int, int >& gap_map1, map<int, int>& gap_map2, ofstream& output, vector<string>& sequences, vector<string>& sequences2, Fasta_vector& file1, Fasta_vector& file2, string& ref_spec,string& ref_spec2, string& outfile, string& files1, string& files2);
int test_num (vector<double>& Correl, double thresh);
double Correlation(vector<double>& sample_a, vector<double>& sample_b);
void inter(vector< vector<double> >& D, vector< vector<double> >& D2, vector< vector<double> >& D_correct, vector< vector<double> >& D_correct2, int& correction, vector<double>& Correl, int simulate, vector<int>& diff1, vector<int>& diff2, vector<string>& sequences, vector<string>& sequences2);
//void print_to_file (vector<double>& Correl, vector< vector<double> >& D, vector< vector<double> >& D_corr, double threshold, vector<string>& sequences, vector<int>& aa_differences, ofstream& output, Fasta_vector& file, string& ref_spec, string& files, int ref);
void print_to_file (vector<double>& Correl, vector< vector<double> >& D, vector< vector<double> >& D_corr, double threshold, vector<string>& sequences, vector<int>& aa_differences, ofstream& output, Fasta_vector& file, string& ref_spec, string& filesi, int ref, double& thresholdR);
void Convert_to_vectors(Fasta_map& file1, Fasta_map& file2, Fasta_vector& vec1, Fasta_vector& vec2);
bool comparison_fabs (double i,double j);
double SD_vf(vector<double>& array, double mean);
int Diff_aa_column(vector<string>& sequences, vector<int>& array);
NUC_AA *allocate_NUC(int size);
int print_splash(string filename);
void Estimate_D(vector<string> seqes, vector< vector<double> >& D, vector< vector<double> >& D_corr, vector<double>& rel_dist, int correction, Blosum& Blos62, vector<int>& aa_differences);
double intra(vector< vector<double> >& D, vector< vector<double> >& D_correct, vector<double>& Correl, int information13, int length_seq, int simulate);
double Li_synonymous(char *seq1, char *seq2, NUC_AA *cero_fold, NUC_AA *two_fold, NUC_AA *four_fold, NUC_AA *six_fold, NUC_AA *purines, NUC_AA *pyrimidines);
TreeTemplate<Node> *create_input_tree(vector<string>& seq_names, vector< string >& sequences);




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  set_hyd
 *  Description:  
 * =====================================================================================
 */
void hydro::set_hyd (){
//	hydro myhydro;

	values['F']=100;
	values['I']=99;
	values['W']=97;
	values['L']=97;
	values['V']=76;
	values['M']=74;
	values['Y']=63;
	values['C']=49;
	values['A']=41;
	values['T']=13;
	values['H']=8;
	values['G']=0;
	values['S']=-5;
	values['Q']=-10;
	values['R']=-14;
	values['K']=-23;
	values['N']=-28;
	values['E']=-31;
	values['P']=-46;
	values['D']=-55;

}		/* -----  end of function set_hyd  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  set_mw
 *  Description:  
 * =====================================================================================
 */
	void MW::set_mw ()
{
	values['F']=165.2;
	values['I']=131.2;
	values['W']=204.2;
	values['L']=131.2;
	values['V']=117.2;
	values['M']=149.2;
	values['Y']=181.2;
	values['C']=121.2;
	values['A']=89.1;
	values['T']=119.1;
	values['H']=155.2;
	values['G']=75.1;
	values['S']=105.1;
	values['Q']=146.2;
	values['R']=174.2;
	values['K']=146.2;
	values['N']=132.1;
	values['E']=147.1;
	values['P']=115.1;
	values['D']=133.1;
}		/* -----  end of function set_mw  ----- */


/*NUC_AA *allocate_NUC(int size){
	NUC_AA *temp;
	int i;

	temp = calloc(1,sizeof(NUC_AA));

	temp->size = size;

	temp->nuc = calloc(size,sizeof(char*));
	temp->amino = calloc(size,sizeof(char*));

	for(i=0;i<size;i++){
		temp->nuc[i] = calloc(4,sizeof(char));
		temp->amino[i] = calloc(6,sizeof(char));
	}

	return temp;
}

void FREE_NUC(int size, NUC_AA *NUC){
	int i;

	for(i=0;i<size;i++){
		if(NUC->nuc[i]!=NULL)
			free(NUC->nuc[i]);
		if(NUC->amino!=NULL)
			free(NUC->amino[i]);
	}	
	if(NUC!=NULL)
		free(NUC);

}
*/
