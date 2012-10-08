#include <iostream>
#include <fstream>
#include <iomanip>

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/SequenceApplicationTools.h>

// From PhylLib:
#include <Phyl/TreeTemplate.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/NonHomogeneousSequenceSimulator.h>
#include <Phyl/HomogeneousSequenceSimulator.h>
#include <Phyl/DRNonHomogeneousTreeLikelihood.h>
#include <Phyl/SubstitutionModelSetTools.h>
#include <Phyl/MarginalAncestralStateReconstruction.h>
#include <Phyl/SequenceSimulationTools.h>
#include <Phyl/SubstitutionModelSetTools.h>
#include <Phyl/Newick.h>
#include<Phyl/models>
#include<Phyl/JCprot.h>
//#include<Phyl/LG08.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/DataTable.h>

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/FileTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/Number.h>

using namespace std;
//#include "fd.h"
using namespace bpp;


std::auto_ptr<DistanceMatrix> ScoreDist(vector<string>& names, vector< string >& sequences);
//int create_seq(TreeTemplate<Node>*& mytree, vector<string>& Tname, vector<string>& tseqs, unsigned int length, const char *mod);
int create_seq(TreeTemplate<Node>*& mytree, vector<string>& Tname, vector<string>& tseqs, unsigned int length, const char *mod, int var, vector<Node *> &nams);
string nuc_to_amino(string nuc);



class codon{

	public:
		map<string, char> convert;
		void initialise();
};
