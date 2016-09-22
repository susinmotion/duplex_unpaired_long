#ifndef LEAFDATA_H
#define LEAFDATA_H
#include <string>
#include <vector>
#include "variant.h"
using namespace std;

class LeafData {
public:
    LeafData(int mCount = 0, int mRevCount=0, string mConsensusFwd="", string mConsensusRev=""): mCount(mCount), mRevCount(mRevCount), mConsensusFwd(mConsensusFwd), mConsensusRev(mConsensusRev), mIsTrash(false){};
    int count();
    int revCount();
    void setCount();
    void setRevCount();
    string consensusFwd();
    string consensusRev();
    vector <Variant*> variants();
    vector <Variant*> superVariants();
    void callConsensus(string sequence, string rev);
    bool isTrash();
    void makeTrash();
    bool hasVariant();
    void setVariants(vector<Variant*> variants);
    void setSuperVariants(vector<Variant*> variants);
    void setGeneLoc(int loc);
    int geneLoc();
private:
    int mGeneLoc;
    int mCount;
    int mRevCount;
    string mConsensusFwd;
    string mConsensusRev;
    vector <Variant*> mVariants;
    vector <Variant*> mSuperVariants;
    bool mIsTrash;
};
#endif
