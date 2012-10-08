#include "BCFasta.h"

int check_alpha(char val, string alpha);
void Fasta_vector::print_to_fasta(const char *filename){


	ofstream file1(filename, ios::app);

		file1 << "\nAmino acid sequences alignment for " << filename << endl;
		file1 << "-----------------------------------------------------------------" << endl;	



	for(int i=0;i<names.size();++i){
		file1 << ">" << names[i] << endl;
		file1 << sequences[i] << endl;
	}
	file1.close();

}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  set_aminos
 *  Description:  
 * =====================================================================================
 */
void set_aminos (map<string, char>& abr)
{

	abr["ALA"]='A';
	abr["ARG"]='R';
	abr["ASN"]='N';
	abr["ASP"]='D';
	abr["CYS"]='C';
	abr["GLN"]='Q';
	abr["GLU"]='E';
	abr["GLY"]='G';
	abr["HIS"]='H';
	abr["ILE"]='I';
	abr["LEU"]='L';
	abr["MET"]='M';
	abr["PHE"]='F';
	abr["PRO"]='P';
	abr["SER"]='S';
	abr["THR"]='T';
	abr["TRP"]='W';
	abr["TYR"]='Y';
	abr["VAL"]='V';
	abr["LYS"]='K';


}		/* -----  end of function set_aminos  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Fasta_vactor::check_for_pdb
 *  Description:  
 * =====================================================================================
 */
void Fasta_vector::check_for_pdb (string& folder, string& filename){


	string tempdb = folder;
	tempdb += "/";
	tempdb += filename;

	vector<string> tok;

	Tokenize(tempdb.c_str(), tok, ".");

	tok[0] += ".pdb";

//	cerr << "tok[0]=" << tok[0] << endl;

	pdb.filename = tok[0];	


	Read_pdb(pdb.filename);
//exit(-1);


}		/* -----  end of function Fasta_vactor::check_for_pdb  ----- */





/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read_pdb
 *  Description:  
 * =====================================================================================
 */
void Fasta_vector::Read_pdb (string& filename){
	ifstream fil;
	char buffer[2048];


	map<string, char> amin;
	set_aminos(amin);
	fil.open(filename.c_str());
	if(fil.is_open()){

	has_pdb=true;

		chain temp_chain;
		amino temp_amino;

		int res = 1;
		char chain;
		bool firsta=true;
		bool firstc=true;
		while(fil.getline(buffer, 2048)){

			if(buffer[0]=='A'&&buffer[1]=='T'&& buffer[2]=='O'&&buffer[3]=='M'){

				//vector<string> tokens;
				//Tokenize(buffer, tokens, " ");
				//cerr << temp << endl;
				string temp;
				temp += buffer;


				string tem(temp, 22, 4);
				atom temp_atom;
				if(firsta==true){
					res=atoi(tem.c_str());
					//++res;
					firsta=false;
				}
				if(res!=atoi(tem.c_str())){
					temp_chain.res_map[res]=temp_chain.aminos.size();
					temp_chain.aminos.push_back(temp_amino);
					temp_amino.atoms.clear();
					res=atoi(tem.c_str());
					//++res;
				//	firsta=true;
				}
				temp_amino.entry=atoi(tem.c_str());
				tem.clear();
				tem.insert(tem.begin(), temp.begin()+12, temp.begin()+16);
				temp_atom.name=tem;
				tem.clear();
				tem.insert(tem.begin(), temp.begin()+17, temp.begin()+20);
				temp_amino.name=tem;
				temp_amino.letter=amin[tem];
				tem.clear();
				if(firstc==true){
					chain=temp[21];
					firstc=false;
				}
				if(temp[21]!=chain){
					pdb.chains.push_back(temp_chain);
					temp_chain.aminos.clear();
					firstc=true;

				}
	
				temp_chain.name=temp[21];
				tem.insert(tem.begin(), temp.begin()+22, temp.begin()+26);

				temp_amino.res=atoi(tem.c_str());


				tem.clear();
				tem.insert(tem.begin(), temp.begin()+30, temp.begin()+38);
				temp_atom.x=atof(tem.c_str());	
				tem.clear();
				tem.insert(tem.begin(), temp.begin()+38, temp.begin()+46);
				temp_atom.y=atof(tem.c_str());	
				tem.clear();
				tem.insert(tem.begin(), temp.begin()+46, temp.begin()+54);
				temp_atom.z=atof(tem.c_str());	
				temp_amino.atoms.push_back(temp_atom);	





				//		tem.insert(temp.begin()+5, temp.begin()+ 	


			}else if(buffer[0]=='S'&& buffer[1]=='O'&&buffer[2]=='U'){
				
				string str, line;
				str += "ORGANISM_COMMON:";
 size_t found;
				line = buffer;
				found = line.find(str);
			if(found!=string::npos){
				vector<string> tok, tok2;
				Tokenize(line, tok, ":");
				Tokenize(tok[1], tok2, ", ;");
				ref_seq=tok2[0];	
	//ref_seq.erase(remove_if(ref_seq.begin(), ref_seq.end(), isspace), ref_seq.end());

				if(ref_seq[0]==' '){
					ref_seq.erase(0,1);

				}
	

			}


			}

		}

		temp_chain.aminos.push_back(temp_amino);
		pdb.chains.push_back(temp_chain);


	}else{
		return;
	}

	pdb.calc_mean_pos();
fil.close();
}		/* -----  end of function Read_pdb  ----- */



void amino::calculate_mean_pos(){


		mean_pos.x=0.0;		
		mean_pos.y=0.0;		
		mean_pos.z=0.0;		
	for(int i=0;i<atoms.size();++i){
		mean_pos.x+=atoms[i].x;		
		mean_pos.y+=atoms[i].y;		
		mean_pos.z+=atoms[i].z;		

	}
		mean_pos.x/=atoms.size();		
		mean_pos.y/=atoms.size();		
		mean_pos.z/=atoms.size();		

}


void structure::calc_mean_pos(){

	for(int i=0;i<chains.size();++i){
		for(int j=0;j<chains[i].aminos.size();++j){
			chains[i].aminos[j].calculate_mean_pos();
		}
	}

}



int Fasta_vector::get_number_of_sequences(){

	return names.size();

}

void Fasta_vector::clear(){
	filename.clear();
	sequences.clear();
	names.clear();
	Tags.clear();

}

void Fasta_map::clear(){
	filename.clear();
	sequences.clear();
	Tags.clear();

}


int Fasta_vector::check_valid_alignment(const char* Alphabet){

	string alpha;

	if(Alphabet=="DNA"){
		alpha = "ACGT-";

	}else if(Alphabet=="AMINO"){
		alpha = "ACDEFGHIKLMNPQRSTVWY-";

	}else if(Alphabet=="AMINOX"){
		alpha = "ACDEFGHIKLMNPQRSTVWYX-";
	}else{
		cerr << "Error: Alphabet other than DNA or AMINO passed to check_valid_alignment" << endl;
		exit(-1);
	}
	int x;

	x=check_length_conservation();
	//cerr << "x=" << x << "\tline=" << __LINE__ << endl;
	if(x==-1){
		return(-1);
	}
	x=check_composition(alpha);
	//cerr << "x=" << x << "\tline=" << __LINE__ << endl;
	if(x==-1){
		return(-1);
	}

	return 0;

}

int Fasta_vector::getAveGaps(){

	int numGaps=0;
	for(int i=0;i<sequences.size();++i){
		for(int j=0;j<sequences[0].size();++j){
			if(sequences[i][j]=='-'){
				++numGaps;
			}		
		}
	}


	return ceil((double)numGaps/(double)sequences.size());

}


int Fasta_vector::check_composition(string alpha){

	int alpha_length = alpha.length();

	for(int i=0;i<sequences.size();++i){
		for(int j=0;j<sequences[i].length();++j){
			//regex rx((const char *)sequences[i][j]);


			if(check_alpha(sequences[i][j], alpha)){	




				//		if(sequences[i].compare(j,1, alpha, 0, alpha_length)!=0){
				cerr << "Warning: The alignment has a non-alphabet character: " << sequences[i][j] << endl;
				cerr << "These will be treated as gaps." << endl;
				return(-1);
			} 
			}

		}

	}


	int Fasta_map::check_composition(string alpha){

		int alpha_length = alpha.length();
		map<string, string>::iterator i;


		for(i=sequences.begin();i!=sequences.end();++i){
			for(int j=0;j<i->second.length();++j){
				//regex rx((const char *)sequences[i][j]);


				if(check_alpha(i->second[j], alpha)){	




					//		if(sequences[i].compare(j,1, alpha, 0, alpha_length)!=0){
				cerr << "Warning: The alignment has a non-alphabet character: " << i->second[j] << endl;
				cerr << "These will be treated as gaps." << endl;
					return(-1);
				} 
				}

			}

		}


		/*		TreeTemplate<Node>* Fasta_Newick_map::create_Bionj(const char* model, const char* alphabet){


				DistanceMatrix *DS;
				Alphabet *alpha1 = new NucleicAlphabet(); 
				Alphabet *alpha2 = new ProteicAlphabet();
				DiscreteDistribution * rdist = new ConstantDistribution(1.);
				if(alphabet=="DNA"){

				}else if(alphabet=="AMINO"){

				const ProteicAlphabet * alpha = dynamic_cast<const ProteicAlphabet *>(alpha2);
				SubstitutionModel * model = new JTT92(alpha);
				VectorSequenceContainer *vsc = new VectorSequenceContainer(alpha2);
				map<string, string>::iterator it;
				for(it=sequences.begin();it!=sequences.end();++it){
				vsc->addSequence(Sequence(it->first, it->second, alpha));
				}
				VectorSiteContainer * sites = new VectorSiteContainer(alpha2);
				for(it=sequences.begin();it!=sequences.end();++it){
				const Sequence *myseq = vsc->getSequence(it->first);
				sites->addSequence(*myseq, true);


				}
				SiteContainerTools::changeGapsToUnknownCharacters(*sites);
				DistanceEstimation MyDS(model, rdist, sites, 1, true);
				delete sites;
				delete vsc;
				delete alphabet;
				delete rdist;
				DS = MyDS.getMatrix();

				}else{
				cerr << "Error: Can only use DNA or AMINO alphabets." << endl; 
				exit(-1);
				}
				delete alpha1;
				delete alpha2;

				}


*/
		TreeTemplate<Node>* Fasta_Newick_map::Read_Newick(const char* filename){



			TreeTemplate<Node> *tree_for = NULL;    
			Newick * NewickReader = new Newick(false); //No comment allowed!
			try{
				tree_for = NewickReader->read(filename);  
			} catch (Exception e){
				cerr << "Error: Couldn't read a Newick tree from file " << filename << endl;

			}
			delete NewickReader;

			return tree_for;




		}


		int check_alpha(char val, string alpha){

			for(int i=0;i<alpha.length();++i){
				if(val==alpha[i]){
					return 0;
				}
			}

			return 1;

		}



		int Fasta_vector::check_length_conservation(){

			int length = sequences[0].length();
			for(int i=1;i<sequences.size();++i){

				if(sequences[i].length()!=length){
					cerr << "Error: Not all of the sequences are the same length, hence not aligned properly" << endl;
					return(-1);
				}
			}
			return(0);
		}

		int Fasta_map::check_length_conservation(){
			map<string, string>::iterator it;

			int length = sequences.begin()->second.length();
			for(it=sequences.begin();it!=sequences.end();++it){

				if(it->second.length()!=length){
					cerr << "Error: Not all of the sequences are the same length, hence not aligned properly" << endl;
					return(-1);
				}
			}
			return(0);
		}


		void Fasta_map::print_to_fasta(const char *filename){

			ofstream file1(filename, ios::app);

			/*				file1 << "\n\t\t*******************************************************************************************\n\t\t";
							file1 << "* CAPS: Co-Evolution Analysis using Protein Sequences                                     *\n\t\t";
							file1 <<  "* Author: Brian E. Caffrey                                                                 *\n\t\t";
							file1 <<  "* Code for Inter-protein co-evolution clustering: David McNally                                  *\n\t\t";
							file1 <<  "* Evolutionary Genetics and Bioinformatics Laboratory                                    *\n\t\t";
							file1 <<  "* Department of Genetics                                                                 *\n\t\t";
							file1 <<  "* Smurfit Institute of Genetics                                                                  *\n\t\t";
							file1 <<  "* University of Dublin, Trinity College                                                          *\n\t\t";
							file1 <<  "* Mathematical Model: Fares and Travers, Genetics (2006)173: 9 - 23                      *\n\t\t";
							file1 <<  "* Conversion to C and addition of multiple allignment functionality: Brian Caffrey (2010) *\n\t\t";
							file1 <<  "*******************************************************************************************\n";
							*/			file1 << "\nAmino acid sequences alignment for " << filename << endl;
			file1 << "-----------------------------------------------------------------" << endl;	

			map<string, string>::iterator it;

			for(it=sequences.begin();it!=sequences.end();++it){
				file1 << ">" << it->first << endl;
				file1 << it->second << endl;
			}


			file1.close();

		}




		int Fasta_map::get_number_of_sequences(){

			return sequences.size();

		}

		void Read_Fasta_vector(const char *fasta_filename, Fasta_vector& fasta){
			/* things to note: file lines shouldn't be longer than 100000 characters 
			   Files should have the regular Fasta format*/

			ifstream file1(fasta_filename); // read the given file

			if(file1.is_open()!=1){
				cerr << "Error: Couldn't open file " << fasta_filename << endl;
				exit(-1);

			}

			char temp[100000];
			string sequence;
			int first=0;
			fasta.filename= fasta_filename;


			file1.getline(temp, 100000);
			if(temp[0]!='>')
				fasta.Tags.push_back(temp);
			while(temp[0]!='>'){

				file1.getline(temp, 100000);
				if(temp[0]!='>')
					fasta.Tags.push_back(temp);

			}
			//temp.erase(0,1);
			fasta.names.push_back(temp);/* push back the first name*/
			//	fasta.names[0].erase(0,1);

			//		fasta.names[0].erase(0,1);

			/*		while(fasta.names[0][0]==' '){
					fasta.names[0].erase(0, 1);
					}
					*/

			int d=1;

			first=1;
			while(file1.getline(temp, 100000)){

				if(temp[0]=='>'){
					//vector<string>::iterator it;
					//fasta.names.find(temp);
					for(int j=0;j<fasta.names.size();++j){
						if(strcmp(temp, fasta.names[j].c_str())==0){
							char tem[10];
							sprintf(tem, "_%d", d);
							strcat(temp, tem);
							++d;
							break;
						}
					}


					fasta.names.push_back(temp);
					//	fasta.names[fasta.names.size()-1].erase(0,1);

					while(fasta.names[fasta.names.size()-1][0]==' '){
						fasta.names[fasta.names.size()-1].erase(0, 1);
					}

					//	if(first!=0){
					fasta.sequences.push_back(sequence);
					sequence.clear();
					//	}
					//	first=1;
				}else{
					sequence += temp;
				}


			}

			for(int i=0;i<fasta.names.size();++i){
				if(fasta.names[i][0]=='>'){
					fasta.names[i].erase(0, 1);
				}
			}

			fasta.sequences.push_back(sequence);

			file1.close();


		}


		void Read_Fasta_map(const char *fasta_filename, Fasta_map& fasta){
			/* things to note: file lines shouldn't be longer than 100000 characters 
			   Fasta files can have junk at the start, it will be written to a vector<string> called
			   tags
			   Files should have the regular Fasta format*/

			ifstream file1(fasta_filename); // read the given file

			if(file1.is_open()!=1){
				cerr << "Error: Couldn't open file " << fasta_filename << endl;
				exit(-1);

			}
			fasta.filename= fasta_filename;
			char temp[100000];
			string sequence, name;
			int first=0;

			file1.getline(temp, 100000);
			if(temp[0]!='>'){
				fasta.Tags.push_back(temp);
			}


			while(temp[0]!='>'){

				file1.getline(temp, 100000);
				if(temp[0]!='>'){
					fasta.Tags.push_back(temp);
					fasta.reference = temp;
				}
			}
			name = temp;
			name.erase(0,1);
			fasta.reference = name;
			while(name[0]==' '){
				name.erase(0, 1);
			}

			while(file1.getline(temp, 100000)){
				//		cerr << "sequence=" << sequence << endl;
				//		cerr << "temp=" << temp << endl;
				if(temp[0]=='>'){


					//		if(first!=0){
					fasta.sequences[name] = sequence;
					//fasta.sequences[name] = sequence;
					sequence.clear();
					//		}

					name = temp;
					name.erase(0, 1);
					while(name[0]==' '){
						name.erase(0, 1);
					}
					first=1;
				}else{
					sequence += temp;
				}


			}
			fasta.sequences[name] = sequence;

			file1.close();


		}

