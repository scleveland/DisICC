#include<iostream>
#include<vector>
#include<fstream>
#include<map>
#include<Phyl/Newick.h>
#include<Seq/Alphabet.h>
#include<Phyl/BioNJ.h>
#include<Phyl/models>
#include<Phyl/DistanceEstimation.h>
#include<Seq/Sequence.h>
#include<Seq/SequenceContainerTools.h>
#include<Seq/SiteContainerTools.h>
#include<Seq/VectorSequenceContainer.h>
#include"file_manip.h"
//#include<boost/regex.hpp>
//#include<boost/regex.hpp>

using namespace std;
using namespace bpp;



class atom{

	public:
		double x;
		double y;
		double z;
		string name;

};


class amino{
	public:
		int entry;
		int res;
		char letter;
		string name;
		vector<atom> atoms;
		atom mean_pos;
		void calculate_mean_pos();

};

class chain{
	public:
		vector<amino> aminos;
		map<int, int> res_map;
		char name;
};


class structure{

	public:
		vector<chain> chains;
		string filename;
		void calc_mean_pos();
};





class Fasta_vector{

	public:
		string filename;
		string reference;
		string ref_seq;
		vector<string> names;
		vector<string> sequences;
		vector<string> Tags;
		bool has_pdb;
		//		string pdb;
		structure pdb;
		void check_for_pdb (string& folder, string& filename);
		int get_number_of_sequences();
		void print_to_fasta(const char *filename);
		int check_valid_alignment(const char *Alphabet);
		void Read_pdb (string& filename);
		void clear();
		int getAveGaps();
	private:
		int check_length_conservation();
		int check_composition(string alpha);
};



class Fasta_map{

	public:
		string filename;
		string reference;
		map<string, string> sequences;
		map<string, string>::iterator sequence_iterator;
		vector<string> Tags;
		int get_number_of_sequences();
		void print_to_fasta(const char *filename);
		int check_valid_alignment(const char *Alphabet);
		void clear();
	private:
		int check_length_conservation();
		int check_composition(string alpha);

};


class Fasta_Newick_map:public Fasta_map{

	public:
		TreeTemplate<Node> *newick;
		TreeTemplate<Node>* Read_Newick(const char* filename);
		TreeTemplate<Node>* create_Bionj(const char* model, const char* alphabet);

};




void set_aminos (map<string, char>& abr);
void Read_Fasta_map(const char *fasta_filename, Fasta_map& fasta);
void Read_Fasta_vector(const char *fasta_filename, Fasta_vector& fasta);
