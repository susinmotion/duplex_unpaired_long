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
string LeafData::superConsensus(){
    return mSuperConsensus;
}

vector <int> LeafData::variants(){
    return mVariants;
}

//this is still a little icky...
void LeafData::callConsensus(string currentSequence, string rev="fwd"){//check substitutions found in this read against those from other reads of this barcode. Non-matches are listed as unconfirmed
    int count;
    string paradigm;
    if (rev=="rev"){
	count=mRevCount;
	paradigm=mConsensusRev;
    }
    else if (rev=="fwd"){
	count=mCount-mRevCount;
	paradigm=mConsensusFwd;
    }
    else{
	cout<<"super!"<<endl;
	count=min(mCount-mRevCount, mRevCount);
	paradigm=mConsensusFwd;
    }
    if (count>0){
	if (paradigm.length()!=currentSequence.length()){//if there's an indel, ignore it and trash
		makeTrash();
		return;
	}
	for (int i=0; i<currentSequence.length(); ++i){//go through existing substitutions. if an item is there but isn't in the new list, mark it as uncertain
		if (currentSequence[i]!=paradigm[i]){
			paradigm.replace(i,1,"N");
		}
	}
    }
    else{
	paradigm=currentSequence;
    }

    if (rev=="rev"){
	mConsensusRev=paradigm;
    }
    else if (rev=="fwd"){
	mConsensusFwd=paradigm;
    }
    else if (rev=="super"){
	mSuperConsensus=paradigm;
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

