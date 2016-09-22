#include "variants.h"
#include "variant.h"
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <map>
#include "leafdata.h"
using namespace std;

string nucleotides = "ACGTN";
vector <Variant*> bowtieCheckVariants(string sequence, string gene, LeafData* data){
            ofstream posfile;
            posfile.open("posfile.txt", ofstream::app);
	char buffer[128];
	string query="bowtie -k 1 "+gene+" --suppress 1,2,3,5,6,7 -c "+sequence+" 2> /dev/null";
	
	string result="";
	FILE* pipe=popen(query.c_str(), "r");
	if (!pipe) {cout<<"no pipe"<<endl; throw runtime_error("popen() failed"); }
	try{
		while (!feof(pipe)){
			if (fgets(buffer,128,pipe)!=NULL)
				result += buffer;
		}
	}
	catch(...){
		pclose(pipe);
		throw;
	}
	pclose(pipe);
	string::size_type pos=0;
	string variant;
	vector<Variant*> variants;
	int locpos;
	int genpos;
	char targ;
	char act;
	bool still_going=true;
	int genomeLoc;
	if (result!=""){
		genomeLoc=atoi(result.substr(0,result.find("\t")).c_str());
		if (data!=NULL){
			data->setGeneLoc(genomeLoc);
			posfile<<genomeLoc<<endl;
		}
		result=result.substr(result.find("\t")+1,result.length()-result.find("\t"));
		if(result != "\n"){
			while ( still_going){
				pos=result.find(",");
				if (pos ==string::npos){
					still_going=false;
					pos=result.length();
				}

			variant=result.substr(0,pos);
			locpos=atoi(variant.substr(0,variant.find(":")).c_str());
			genpos=locpos+genomeLoc;
			targ=variant[variant.find(">")-1];
			act=variant[variant.find(">")+1];
			Variant* v= new Variant(genpos,act,targ);
			variants.push_back(v);
			result.erase(0, pos+ 1);
			}
        	}
	}
	posfile.close();
return variants;
}

	
