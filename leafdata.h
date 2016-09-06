#ifndef LEAFDATA_H
#define LEAFDATA_H
#include <string>
#include <vector>
using namespace std;

class LeafData {
public:
    LeafData(int mCount = 0, int mRevCount=0, string mConsensusFwd="", string mConsensusRev=""): mCount(mCount), mRevCount(mRevCount), mConsensusFwd(mConsensusFwd), mConsensusRev(mConsensusRev){};
    int count();
    int revCount();
    void setCount();
    void setRevCount();
    string consensusFwd();
    string consensusRev();
    vector <int> variants();
    void callConsensus(string sequence, bool rev);
    bool isTrash();
    void makeTrash();
    bool hasVariant();
    void setVariants(vector<int> variants);

private:
    int mCount;
    int mRevCount;
    string mConsensusFwd;
    string mConsensusRev;
    vector <int> mVariants;
    bool mIsTrash;
};
#endif
