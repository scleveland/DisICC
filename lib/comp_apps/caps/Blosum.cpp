/*
 * =====================================================================================
 *
 *       Filename:  Blosum.cpp
 *
 *    Description:  File which has the Blosum matrices hard coded
 *
 *        Version:  1.0
 *        Created:  17/09/2010 11:55:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Brian E. Caffrey (), brian.caffrey@gmail.com
 *        Company:  Smurfit Institute of Genetics
 *
 * =====================================================================================
 */

#include"Blosum.h"

Blosum Blosum62(){

	Blosum temp62;

	vector<int> temp_val;

	string str= "A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y  -";
	vector<string> tokens;
	Tokenize(str, tokens, " ");

	for(int i=0;i<tokens.size();++i){
		temp62.index[tokens[i]]=i;
	}

	int myints[] = {4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2,0};
	temp_val.assign(myints, myints+21);
	temp62.values.push_back(temp_val);

	int myints1[] = {0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2,0};
	temp_val.assign(myints1, myints1+21);
	temp62.values.push_back(temp_val);
	int myints2[] = {-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3,0};
	temp_val.assign(myints2, myints2+21);
	temp62.values.push_back(temp_val);
	int myints3[] = {-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2,0};
	temp_val.assign(myints3, myints3+21);
	temp62.values.push_back(temp_val);
	int myints4[] = {-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3,0};
	temp_val.assign(myints4, myints4+21);
	temp62.values.push_back(temp_val);
	int myints5[] = {0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3,0};
	temp_val.assign(myints5, myints5+21);
	temp62.values.push_back(temp_val);
	int myints6[] = {-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2,0};
	temp_val.assign(myints6, myints6+21);
	temp62.values.push_back(temp_val);
	int myints7[] = {-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1,0};
	temp_val.assign(myints7, myints7+21);
	temp62.values.push_back(temp_val);
	int myints8[] = {-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2,0};
	temp_val.assign(myints8, myints8+21);
	temp62.values.push_back(temp_val);
	int myints9[] = {-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1,0};
	temp_val.assign(myints9, myints9+21);
	temp62.values.push_back(temp_val);
	int myints10[] = {-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1,0};
	temp_val.assign(myints10, myints10+21);
	temp62.values.push_back(temp_val);
	int myints11[] = {-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2,0};
	temp_val.assign(myints11, myints11+21);
	temp62.values.push_back(temp_val);
	int myints12[] = {-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3,0};
	temp_val.assign(myints12, myints12+21);
	temp62.values.push_back(temp_val);
	int myints13[] = {-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1,0};
	temp_val.assign(myints13, myints13+21);
	temp62.values.push_back(temp_val);
	int myints14[] = {-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2,0};
	temp_val.assign(myints14, myints14+21);
	temp62.values.push_back(temp_val);
	int myints15[] = {1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2,0};
	temp_val.assign(myints15, myints15+21);
	temp62.values.push_back(temp_val);
	int myints16[] = {0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2,0};
	temp_val.assign(myints16, myints16+21);
	temp62.values.push_back(temp_val);
	int myints17[] = {0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1,0};
	temp_val.assign(myints17, myints17+21);
	temp62.values.push_back(temp_val);
	int myints18[] = {-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2,0};
	temp_val.assign(myints18, myints18+21);
	temp62.values.push_back(temp_val);
	int myints19[] = {-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7,0};
	temp_val.assign(myints19, myints19+21);
	temp62.values.push_back(temp_val);
	int myints20[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	temp_val.assign(myints20, myints20+21);
	temp62.values.push_back(temp_val);

	return temp62;

}

void Blosum::print(const char * filename){

ofstream file1(filename);


if(file1.is_open()){
}else{
	cerr << "Error: couldn't open file " << filename << endl;
	exit(-1);
}


	map<string, int>::iterator myit, myit1;
	for(myit=index.begin();myit!=index.end();++myit){
		file1 << myit->first << "\t";		
	}
		file1 << endl;
	for(myit=index.begin();myit!=index.end();++myit){
	for(myit1=index.begin();myit1!=index.end();++myit1){
		file1 << values[index[myit->first]][index[myit1->first]] << "\t";
	}
	file1 << endl;
	}	

file1.close();

}

