#include "trie.h"
#include <string>
#include <map>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <iostream>
using namespace std;

string reverseComplement(string sequence){
    map <char, char> complements;
    complements['A']='T';
    complements ['T']='A';
    complements ['G']='C';
    complements ['C']='G';
    complements ['N']='N';

    string reverseComplement="";
 
    for (int i = sequence.length()-1; i>=0; --i){
        reverseComplement = reverseComplement+complements[sequence[i]];
    }
    return reverseComplement;
}

map <string, vector <string> > readConfig(string filename){ //read config file into map of user defined variables
    map <string, vector<string> > userDefinedVariables;
    userDefinedVariables.clear();
    ifstream infile(filename.c_str());
    if (infile.is_open()){
        string key;
        string allValues;
        string value;
        while (infile >> key >> allValues ){
            vector <string> values;
            cout <<key<<" "<<allValues<<endl;
            stringstream ss(allValues);
            while (getline (ss, value, ',')){
                values.push_back(value);
            }

            pair <string, vector<string> > varPair = make_pair (key, values);
            userDefinedVariables.insert(varPair);
        }
    }
    else {
        cout<<"error opening config file "<<filename<<endl;
    }
    return userDefinedVariables;
}
    
Trie* readFileIntoTrie(string filename, int no_of_files, char* pipes[]){//set constants based on config file
    map <string, vector <string> > userDefinedVariables=readConfig(filename);
    cout<<" done reading config "<<endl;
    const int BARCODE_LENGTH =atoi(userDefinedVariables["BARCODE_LENGTH"][0].c_str() );
    const vector <string> GENES = userDefinedVariables["GENES"];
    vector <string> FORWARD_ALIGN_SEQ=userDefinedVariables["FORWARD_ALIGN_SEQ"];
    vector <string> REVERSE_ALIGN_SEQ=userDefinedVariables["REVERSE_ALIGN_SEQ"];
    vector <string> TARGET=userDefinedVariables["TARGET"];
    vector <string> PHASE_SHIFTS=userDefinedVariables["PHASE_SHIFTS_REV_TO_FORWARD"];
    const vector <string> FILENAMES =userDefinedVariables["FILENAMES"];
    vector <int> TARGET_LENGTHS;
    vector <int> THRESHOLDS_OF_IMPORTANCE;

    for (int i=0; i<userDefinedVariables["THRESHOLD_OF_IMPORTANCE"].size(); ++i){
        THRESHOLDS_OF_IMPORTANCE.push_back(atoi(userDefinedVariables["THRESHOLD_OF_IMPORTANCE"][i].c_str()));
    }

    for (int i=0; i<TARGET.size(); ++i){
        TARGET_LENGTHS.push_back(TARGET[i].length());
    }
    
    
    int numberOfROIs=FORWARD_ALIGN_SEQ.size();
    int numberOfPhases=atoi(userDefinedVariables["MAX_PHASE"][0].c_str() )+1;

    map<int, int> empty_map;
    vector <map <int, int> > PHASE_MAPS=vector <map <int, int> >(numberOfROIs, empty_map);
    string keyvalue;
    for (int i=0; i<PHASE_SHIFTS.size(); ++i){
        stringstream ss(PHASE_SHIFTS[i]);
        while (getline (ss, keyvalue, '|')){
            int key =atoi((keyvalue.substr(0, keyvalue.find(":"))).c_str());
            int value=atoi((keyvalue.substr(keyvalue.find(":")+1,keyvalue.length())).c_str());
            PHASE_MAPS[i].insert(make_pair(key, value));
        }

    }
    Trie* trie = new Trie;
    trie->setThresholdROIPhaseGenesBarcodelenTargetlen( THRESHOLDS_OF_IMPORTANCE, numberOfROIs, numberOfPhases, GENES, BARCODE_LENGTH, TARGET_LENGTHS);
    
int count2=0;
    int count=0;
    string sequence;
    string ROI;
    string barcode;
    string throwoutstring;
    string target="";
    int phase=0;
    int ROINumber=0;
    int indexForwardAlign;
    int indexReverseAlign;
    ifstream file1, file2;
    for (int i=1; i<no_of_files; i+=2){
	cout<<pipes[i]<<" "<<pipes[i+1]<<endl;
	file1.open(pipes[i], ifstream::in);
    	file2.open(pipes[i+1], ifstream::in);
	while (getline(file1,throwoutstring)){//read sequence. 4 lines is a read. 2nd line has sequence
            count++;
	    if (count%100000==0){
		cout<<count<<endl;
	    }
            file1>>sequence;
	    if (sequence.substr(BARCODE_LENGTH,5)==FORWARD_ALIGN_SEQ[0]){
		barcode=sequence.substr(0,BARCODE_LENGTH);
		ROI=sequence.substr(BARCODE_LENGTH+5, sequence.length()-BARCODE_LENGTH-4);
		trie->addBarcode(ROINumber, phase,barcode,ROI,target);
		}
	    getline(file2,throwoutstring);
	    file2>>sequence;
	    if (sequence.substr(BARCODE_LENGTH,5)==REVERSE_ALIGN_SEQ[0]){       
                barcode=sequence.substr(0,BARCODE_LENGTH);
                ROI=sequence.substr(BARCODE_LENGTH+5, sequence.length()-BARCODE_LENGTH-4);
		trie->addBarcode(ROINumber, phase,barcode,ROI, target, "rev");                   

	     }
	    
	    getline(file1,throwoutstring);
            getline(file1,throwoutstring);
            getline(file1,throwoutstring);                 
            getline(file2,throwoutstring);
            getline(file2,throwoutstring);
            getline(file2,throwoutstring);
       }
//  }
// else {
//    cout<<"Error opening file "<< userDefinedVariables["FILENAME"][i]<<endl;
//    }
    }
  cout<<"made the trie"<<endl;
   return trie;
}

