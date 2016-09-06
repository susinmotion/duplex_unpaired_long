#include "variants.h"
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
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
