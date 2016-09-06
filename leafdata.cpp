#include <vector>
#include <stack>
#include <algorithm>
#include "leafdata.h"
#include <iostream>
using namespace std;

int LeafData::count() {
    return mCount;
}
int LeafData::revCount() {
    return mRevCount;
}

void LeafData::setCount() {
    mCount++; 
}
void LeafData::setRevCount(){
   mRevCount++;
}

string LeafData::consensusFwd(){
    return mConsensusFwd;
}

string LeafData::consensusRev(){
    return mConsensusRev;
}
vector <int> LeafData::variants(){
    return mVariants;
}

//this is still a little icky...
void LeafData::callConsensus(string currentSequence, bool rev=false){//check substitutions found in this read against those from other reads of this barcode. Non-matches are listed as unconfirmed

	
    if (rev==true){
	if (mRevCount>0){
	    if (mConsensusRev.length()!=currentSequence.length()){//if there's an indel, ignore it and trash
		makeTrash();
	    	return;
	    }	
            for (int i=0; i<currentSequence.length(); ++i){//go through existing substitutions. if an item is there but isn't in the new list, mark it as uncertain
	    	if (currentSequence[i]!=mConsensusRev[i]){
	    	    mConsensusRev.replace(i, 1,"N");
	    	}
	    }
	}
	else{
	    mConsensusRev=currentSequence;
	}
    }
    else{
	if (mCount-mRevCount>0){
	    if (mConsensusFwd.length()!=currentSequence.length()){//if there's an indel, ignore it and trash
                makeTrash();
                return;
            }               
	    for (int i=0; i<currentSequence.length(); ++i){//go through existing substitutions. if an item is there but isn't in the new list, mark it as uncertain
                if (currentSequence[i]!=mConsensusFwd[i]){
                    mConsensusFwd.replace(i, 1,"N");
                }
            }
	}
    	else {//if this is the first read, no check is necessary
            mConsensusFwd=currentSequence;
	}
    }
}
void LeafData::setVariants(vector<int> variants){
    mVariants=variants;
}

bool LeafData::isTrash(){
    return mIsTrash;
}

void LeafData::makeTrash(){
    mIsTrash=true;
}

