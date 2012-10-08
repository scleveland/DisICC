/*
 * =====================================================================================
 *
 *       Filename:  Blosum.h
 *
 *    Description:  header file for Blosum.cpp
 *
 *        Version:  1.0
 *        Created:  21/09/2010 14:01:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Brian E. Caffrey (), brian.caffrey@gmail.com
 *        Company:  Smurfit Institute of Genetics
 *
 * =====================================================================================
 */

#include"file_manip.h"
#include<iostream>
#include<vector>
#include<map>

using namespace std;

class Blosum{

	public:
		vector< vector<int> > values;
		map<string, int> index;
		void print(const char *filename);
	private:


};

Blosum Blosum62();
