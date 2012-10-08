
/*

Copyright 2010 Brian Caffrey, Tom Williams, Mario Fares.


this file is part of Clusterfunc.

    Clusterfunc is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Clusterfunc is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Clusterfunc.  If not, see <http://www.gnu.org/licenses/>.


*/



#include<stdio.h>
#include<cstdlib>
#include<iostream>
#include<string.h>
#include<fstream>
#include<sstream>
#include<dirent.h>//readdir etc
#include<sys/stat.h>
#include<unistd.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include<map>
#include<vector>



using namespace std;


void Pwd(string& pwd);
void Tokenize(const string& str, vector<string>& tokens, const string& delimiters);
int tab_delim_to_map(int col, const char *filename, map<string, vector< vector<string> > >& mymap);
int comma_sep_to_map(int col, const char *filename, map<string, vector< vector<string> > >& mymap);
vector<string> Folder_to_vector(const char *dir);


