#include "variants.h"
#include "variant.h"
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <map>

using namespace std;

string nucleotides = "ACGTN";

int hashVariants (string trio, char shift){//hash pos/nucleotide to int
    int variantHash=nucleotides.find(trio[0])*125+nucleotides.find(trio[1])*25+nucleotides.find(trio[2])*5+nucleotides.find(shift);
    return variantHash;
}

pair<string,char> unhashVariants (int variantHash){//unhash int into pos/nucleotide
    string trio="";
    trio+=nucleotides[variantHash/125];
    trio+=nucleotides[(variantHash%125)/25];
    trio+=nucleotides[(variantHash%25)/5];
    char shift=nucleotides[variantHash%5];
    return pair<string, char>(trio, shift);
}

void checkVariants(LeafData* pCurrentData){
    string currentFwd=pCurrentData->consensusFwd();
    string currentRev=pCurrentData->consensusRev();
    vector <int> currentVariants;
    if (currentRev.length()!=currentFwd.length()){
	pCurrentData->makeTrash();
    }
    else{
	for (int i=1; i<currentRev.length()-1; ++i){
	    if (currentFwd[i]!=currentRev[i]){
    currentVariants.push_back(hashVariants(currentFwd.substr(i-1,3),currentRev[i]));
	    }
	}
    }
    pCurrentData->setVariants(currentVariants);
}

vector<Variant*> bowtieCheckVariants(string sequence, string gene){
    char buffer[128];
    string query="bowtie -k 1 "+gene+" --suppress 1,3,5,6,7 --best -c "+sequence+" 2> /dev/null";
    FILE* pipe=popen(query.c_str(), "r");
    string result="";
    if (!pipe) {throw runtime_error("popen() failed"); }
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
    string subpos;
    int spos;
    char targ;
    char act;
    string strand;
    int alignpos;
    if (result != "\n"&&result!=""){
        istringstream iss(result);
        iss>>strand>>alignpos;
        string rest;
        iss>>rest;
        if (rest!=""){
                while ( 1 ){
                    subpos=rest.substr(0,rest.find(":"));
                    istringstream(subpos)>>spos;
                    targ=rest[rest.find(">")-1];
                    act=rest[rest.find(">")+1];
                    if( strand=="+"){
                        spos=alignpos+spos;
                    }
                    else {
                        spos=alignpos+sequence.length()-spos-1;
                    }
                    Variant* v= new Variant(spos,act,targ);
                    variants.push_back(v);
                    pos=rest.find(",");
                    if (pos==string::npos){
                        break;
                    }
                    else {
                        rest.erase(0, pos+1);
                    }

                }
        }
    }
	return variants;
}

