#include"create.h"
#include<Phyl/DistanceEstimation.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  create_seq
 *  Description:  create sequences for simultaionts 
 * =====================================================================================
 */
int create_seq(TreeTemplate<Node>*& mytree, vector<string>& Tname, vector<string>& tseqs, unsigned int length, const char *mod, int var, vector<Node *> &nams)
{


	DiscreteDistribution  *rDist = new ConstantDistribution(1.);
	Alphabet * alphabet = new ProteicAlphabet(); 
	//Alphabet * alphabet = new ProteicAlphabet();
	vector<Node *> nodes; 

	unsigned int simlength = 0;


	if(var==0){
		nodes = mytree->getNodes();
	}else{
		nodes = nams;
	}
	SubstitutionModel  * model =NULL;
	//		SubstitutionModel * model = new JCprot(dynamic_cast<const ProteicAlphabet *>(alphabet)); 
	if(strcmp(mod, "JC")==0){
		const ProteicAlphabet * alpha = dynamic_cast<const ProteicAlphabet *>(alphabet);
ProteinFrequenciesSet *fSet = new FixedProteinFrequenciesSet(dynamic_cast<ProteicAlphabet *>(alphabet));

		SubstitutionModel * model1 = new JCprot(dynamic_cast<ProteicAlphabet *>(alphabet), fSet);
//		model = new JCProt(dynamic_cast<ProteicAlphabet *>(alpha)); 
		//cerr << "line=" << __LINE__ << endl;
	model = model1->clone();
		//model = new JCProt(alpha); 
	//	model = JCProt(alpha); 
		delete model1;
		simlength = length;
	}else if(strcmp(mod, "JTT")==0){
		model=NULL;
		model = new JTT92(dynamic_cast<ProteicAlphabet *>(alphabet)); 
		//model = new LG08(dynamic_cast<ProteicAlphabet *>(alphabet)); 
		simlength = length;
	}else{
		cerr << "WTF" << endl;
		exit(-1);
	}
	const Alphabet *al = model->getAlphabet();
	Vdouble test = model->getFrequencies();

	FullFrequenciesSet* fSet = new FullFrequenciesSet(model->getAlphabet());
	fSet->setNamespace("anc.");
	fSet->setFrequencies(model->getFrequencies());
	SubstitutionModelSet *modelSet;
	modelSet = SubstitutionModelSetTools::createHomogeneousModelSet(dynamic_cast<SubstitutionModel*>(model->clone()), fSet, mytree);

	NonHomogeneousSequenceSimulator seqsim(modelSet, rDist, mytree);
	//	HomogeneousSequenceSimulator seqsim(model, rDist, mytree);
	//	NonHomogeneousSequenceSimulator *seqsim = new NonHomogeneousSequenceSimulator(model, rDist, mytree);
	std::auto_ptr<SiteContainer> sites(seqsim.simulate(simlength));

	Tname = sites->getSequencesNames();
	for(unsigned int i=0;i<Tname.size();++i){
		const Sequence myseq=sites->getSequence(Tname[i]);
		if(strcmp(mod, "JC")==0){
			tseqs.push_back(myseq.toString());
		}else{
			tseqs.push_back(myseq.toString());
		}
		//delete myseq;
	}

	delete modelSet;
	//delete sites;
	//delete seqsim;
	delete model;
	delete alphabet;
	//delete trantree;
	delete rDist;
	return (0);
}

/* -----  end of function create_seq  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ScoreDist
 *  Description:  Get the Distance matrix for a set of alignments according to JTT model 
 * =====================================================================================
 */

std::auto_ptr<DistanceMatrix> ScoreDist(vector<string>& names, vector< string >& sequences){

	FILE *stream;
	stream = freopen("/dev/null", "w", stdout);
	//DistanceMatrix *DS;
	std::auto_ptr<DistanceMatrix> DS;
	Alphabet *alphabet = new ProteicAlphabet();
	//DiscreteDistribution * rdist = new ConstantDistribution(1.);
	DiscreteDistribution * rdist = new GammaDiscreteDistribution(4, 0.5);
	///ConstantDistribution rdist(1.);
	const ProteicAlphabet * alpha = dynamic_cast<const ProteicAlphabet *>(alphabet);
	SubstitutionModel * model = new JTT92(alpha);
	VectorSequenceContainer *vsc = new VectorSequenceContainer(alpha);

	//SequenceContainer *sc;
	//sc->createEmptyContainer();


	int nam_siz=names.size();	
	for(unsigned int i=0;i<nam_siz;++i){
		std::string mystr(sequences[i]);
		std::string tem(names[i]);
		vsc->addSequence(Sequence(names[i], sequences[i], alphabet));
		//		sc->addSequence(Sequence(names[i], sequences[i], alphabet), false);
	}

	VectorSiteContainer * sites = new VectorSiteContainer(alphabet);

	for(unsigned int i=0;i<nam_siz;++i){
		const Sequence myseq=  vsc->getSequence(names[i]);
		sites->addSequence(myseq, true);	
	}

	model->setFreqFromData(*sites);


	SiteContainerTools::changeGapsToUnknownCharacters(*sites);
	DistanceEstimation  MyDS(model, rdist, sites, 1, true);

	delete sites;
	delete vsc;
	delete model;
	DS = (std::auto_ptr<DistanceMatrix>)(MyDS.getMatrix());
	delete rdist;
	//delete MyDS;
	//delete alpha;
	delete alphabet;

	return DS;
}

/* -----  end of function ScoreDist  ----- */

string nuc_to_amino(string nuc){


	string prot;
	codon mycod;
	mycod.initialise();

	for(unsigned int i=0;i<nuc.size();i+=3){

		string temp;
		temp += nuc[i];
		temp += nuc[i+1];
		temp += nuc[i+2];

		prot += mycod.convert[temp];

	}


	return prot;
}



void codon::initialise(){


	convert["ATG"]='M';
	convert["TAA"]='-';
	convert["TAG"]='-';
	convert["TGA"]='-';
	convert["TGG"]='W';

	convert["TTT"]='F';	
	convert["TTC"]='F';	
	convert["TAT"]='Y';	
	convert["TAC"]='Y';	
	convert["CAT"]='H';	
	convert["CAC"]='H';	
	convert["CAA"]='Q';	
	convert["CAG"]='Q';	
	convert["AAA"]='K';	
	convert["AAG"]='K';	
	convert["GAT"]='D';	
	convert["GAC"]='D';	
	convert["GAA"]='E';	
	convert["GAG"]='E';	
	convert["AAT"]='N';	
	convert["AAC"]='N';	
	convert["TGT"]='C';	
	convert["TGC"]='C';	

	convert["ATT"]='I';	
	convert["ATC"]='I';	
	convert["ATA"]='I';	
	convert["GTA"]='V';	
	convert["GTT"]='V';	
	convert["GTC"]='V';	
	convert["GTG"]='V';	
	convert["CCA"]='P';	
	convert["CCT"]='P';	
	convert["CCC"]='P';	
	convert["CCG"]='P';	
	convert["ACA"]='T';	
	convert["ACT"]='T';	
	convert["ACC"]='T';	
	convert["ACG"]='T';	
	convert["GCA"]='A';	
	convert["GCT"]='A';	
	convert["GCC"]='A';	
	convert["GCG"]='A';	
	convert["GGA"]='G';	
	convert["GGT"]='G';	
	convert["GGC"]='G';	
	convert["GGG"]='G';	

	convert["TTA"]='L';	
	convert["TTG"]='L';	
	convert["CTA"]='L';	
	convert["CTT"]='L';	
	convert["CTC"]='L';	
	convert["CTG"]='L';	
	convert["TCA"]='S';	
	convert["TCT"]='S';	
	convert["TCC"]='S';	
	convert["TCG"]='S';	
	convert["AGT"]='S';	
	convert["AGC"]='S';	
	convert["AGA"]='R';	
	convert["AGG"]='R';	
	convert["CGA"]='R';	
	convert["CGT"]='R';	
	convert["CGC"]='R';	
	convert["CGG"]='R';	


}

