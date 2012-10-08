

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

#include"file_manip.h"

vector<string> Folder_to_vector(const char *dir){
	DIR *dp;
	struct dirent *entry;
	struct stat statbuf;
	string pwd;

	Pwd(pwd);

	vector<string> files;

	if ((dp = opendir(dir)) == NULL) {
		fprintf(stderr, "Cannot open directory: %s\n", dir);
		exit(-1);
	}
	chdir(dir);
	while ((entry = readdir(dp)) != NULL) {
		lstat(entry->d_name, &statbuf);

		if (S_ISREG(statbuf.st_mode)) {
			if(*entry->d_name!='.'){
				files.push_back(entry->d_name);
			}
		} else if (S_ISDIR(statbuf.st_mode)) {
			/* Maybe traverse inside see below on how to do. */
		}

	}
	chdir(pwd.c_str());
	closedir(dp);
	return files;
}

void Pwd(string& pwd){
	char temp[20000];
	getcwd(temp, 20000);
	strcat(temp, "/");

	pwd = temp;

}      



/* calculate filesize */
/*int calc_bufsize(char *filename)
{
        struct stat st;

        stat(filename, &st);
        return ((int)st.st_size);
}
*/

int tab_delim_to_map(int col, const char *filename, map<string, vector< vector<string> > >& mymap){

        fstream file(filename);	

	if(file.is_open()==false){
		cerr << "Error: couldn't open file " << filename << endl;
	}

	char buffer[100000];	

	while(file.getline(buffer, 100000)){

		vector<string> temp_vec;
		Tokenize(buffer, temp_vec, "\t");
		if(mymap.find(temp_vec[col])==mymap.end()){
			vector< vector<string> > temptemp;
			temptemp.push_back(temp_vec);
		mymap[temp_vec[col]] = temptemp;
		}else{
		mymap[temp_vec[col]].push_back(temp_vec);
		}
	}

/*map<string, vector< vector<string> > >::iterator ip;
for(ip=mymap.begin();ip!=mymap.end();++ip){
	for(int i=0;i<ip->second.size();++i){
		for(int j=0;j<ip->second[i].size();++j){
		cout << ip->first << "\t" << ip->second[i][j] << endl;
	}

	}
}
*/
	file.close();
}


int comma_sep_to_map(int col, const char *filename, map<string, vector< vector<string> > >& mymap){

        fstream file(filename);	
	char buffer[100000];	
	if(file.is_open()==false){
		cerr << "Error: couldn't open file " << filename << endl;
	}
	while(file.getline(buffer, 100000)){
		vector<string> temp_vec;
		Tokenize(buffer, temp_vec, ",");
		if(mymap.find(temp_vec[col])==mymap.end()){
			vector< vector<string> > temptemp;
			temptemp.push_back(temp_vec);
		mymap[temp_vec[col]] = temptemp;
		}else{
		mymap[temp_vec[col]].push_back(temp_vec);
		}


	}
	file.close();

}

int grep(const char* filename, const char* pattern)
{

        fstream fp;
	fp.open(filename);	

                if(fp.is_open()){}else
                {
                        printf("Failed to open file: %s\n", filename);
                        return -2;
                }

               // int BUFSIZE = calc_bufsize(filename);

                /* read ENTIRE file into buf[] */
               // char buf[BUFSIZE];
               // fread(&buf, sizeof(char), BUFSIZE, fp);

                /* search buf for word (case sensitive) */

		char temp[1000000];
int line=1;
		while(fp.getline(temp, 1000000)){


                char *ans = strstr(temp, pattern);

                /* word found, print line */
                if(ans != NULL)
                        printf("%d: %s\n", line, temp);

		++line;

		}

                /* word not found, do nothing */

                fp.close();
        return 0;
}

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters){
	//        const string& delimiters = " \t";
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0); 
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	while (string::npos != pos || string::npos != lastPos)
	{   
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}   
}

