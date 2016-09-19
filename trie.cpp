#include "trie.h"
#include "node.h"
#include "variants.h"
#include <iostream>
#include <vector>
#include <string>
#include <stack>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <typeinfo>
#include <iomanip>
#include <set>
using namespace std;

Node* Trie::pRootPointer(){
    return mRootPointer;
}

void Trie::addBarcode(int ROINumber, int phase, string barcode, string sequence, string target, bool rev){
    Node* pCurrentNode = mRootPointer;
    if (barcode.find('N') != -1){ return;}
    if ( barcode.length() == 0 ){
        pCurrentNode->initializeLeafData(mNumberOfROIs, mNumberOfPhases); // an empty barcode
        return;
    }
    for ( int i = 0; i < barcode.length(); i++ ){//go through barcode base by base. if it's in the trie, continue. if not, add it. 
        Node* pChildNode = pCurrentNode->findChild(barcode[i]);
        if ( pChildNode != NULL ){
        pCurrentNode = pChildNode;
	}
        else {
            Node* pTmp = new Node();
            pTmp->setContent(barcode[i]);
	    pCurrentNode->appendChild(pTmp);
            pCurrentNode = pTmp;
        }
        if ( i == barcode.length() - 1 ){//if we are at the end of the barcode, check for variants.
	    if (pCurrentNode->leafData().empty()){
                pCurrentNode->initializeLeafData(mNumberOfROIs, mNumberOfPhases);
            }
            LeafData* pCurrentData= pCurrentNode->leafData()[ROINumber][phase];
            if (pCurrentData==NULL){
                pCurrentData=new LeafData;
            }
	    pCurrentData->callConsensus(sequence, rev);
            pCurrentData->setCount();
	    if (rev==true){
		pCurrentData->setRevCount();
	    }
            pCurrentNode->setLeafData(ROINumber, phase, pCurrentData);
            if (pCurrentNode->leafData()[ROINumber][phase]->count()==mThresholdsOfImportance[0]){//if there are enough reads, add pointer to list of important nodes for output later
                addImportantNode(pCurrentNode, ROINumber, phase);
            }
            if (pCurrentNode->leafData()[ROINumber][phase]->count()>=mThresholdsOfImportance[0] && pCurrentNode->leafData()[ROINumber][phase]->count()<=200){
                mCounts[ROINumber][phase][pCurrentNode->leafData()[ROINumber][phase]->count()]++;
            }
        }
    }
}

void Trie::printCounts(){
    ofstream outfile2;
    outfile2.open("total_counts.txt");
    int counts [mNumberOfROIs];
    for (int i=0; i<mNumberOfROIs; ++i){
        ofstream outfile;
        string outfilename="counts_"+mGenes[i]+".txt";
        outfile.open(outfilename.c_str());
        counts[i]=0;
        outfile<<"ROI "<<mGenes[i]<<endl;
        outfile2<<"ROI "<<mGenes[i]<<" ";
        cout<<"ROI "<<mGenes[i]<<endl;
        for (int j=0; j<mNumberOfPhases; ++j){
            outfile<<" phase "<<j<<":"<<endl;
            for (int k=0; k<mCounts[i][j].size()-1; ++k){
	        if (mCounts[i][j][k]-mCounts[i][j][k+1]>0){
                    outfile<<"  "<<mCounts[i][j][k]-mCounts[i][j][k+1]<<" nodes were found exactly "<<k<<" times"<<endl;
                    counts[i]+=mCounts[i][j][k]-mCounts[i][j][k+1];
		}
            }
        }
        outfile2<<counts[i]<<endl;
    }
}
void Trie::setThresholdROIPhaseGenesBarcodelenTargetlen(vector <int> threshold, int numberOfROIs, int numberOfPhases, vector<string>genes, int barcodeLength, vector <int> targetLength){
    mThresholdsOfImportance=threshold;
    mNumberOfROIs= numberOfROIs;
    mNumberOfPhases = numberOfPhases;
    mGenes = genes;
    mTargetLength= targetLength;
    set <Node*> empty_set;
    mAs=vector<vector<int> > (mNumberOfROIs, vector <int>(mNumberOfPhases,0));
    mGs=vector<vector<int> > (mNumberOfROIs, vector <int>(mNumberOfPhases,0));
    mCs=vector<vector<int> > (mNumberOfROIs, vector <int>(mNumberOfPhases,0));
    mTs=vector<vector<int> > (mNumberOfROIs, vector <int>(mNumberOfPhases,0));
    mNs=vector<vector<int> > (mNumberOfROIs, vector <int>(mNumberOfPhases,0));
     mImportantNodes=vector <vector <set <Node*> > >(mNumberOfROIs, vector<set<Node* > >(mNumberOfPhases, empty_set));
    mBarcodeLength=barcodeLength;
    mCounts=vector< vector <vector<int> > >(mNumberOfROIs, vector< vector <int> >(mNumberOfPhases,vector <int>(200,0)));
}

vector< vector< set <Node*> > >Trie::importantNodes(){
    return mImportantNodes;
}

void Trie::addImportantNode(Node* pImportantNode, int ROINumber, int phase){
    mImportantNodes[ROINumber][phase].insert(pImportantNode);
}

void Trie::populateVariants(int threshold){
    mVariantsCount= vector <vector<int> >(mNumberOfROIs, vector<int>(mNumberOfPhases,  0));
    mNodesChecked= vector <vector<int> >(mNumberOfROIs, vector<int>(mNumberOfPhases,  0));
    mVariants=vector <vector <vector<int> > >(mNumberOfROIs, vector<vector< int> >(mNumberOfPhases, vector<int>(625, 0)));
              
    cout<<mNumberOfROIs<< " rois"<<endl;
    cout<<mNumberOfPhases<<" phases"<<endl;
    int totalImportantNodes=0;
    for (int i=0; i<mNumberOfROIs; ++i){
	for (int j=0; j<mNumberOfPhases; ++j){
	    set <Node*> currentImportantNodes=mImportantNodes[i][j];
	    cout<<currentImportantNodes.size()<<" important nodes"<<endl;
	    for (set <Node*>::iterator it=currentImportantNodes.begin(); it!=currentImportantNodes.end(); ++it){
		LeafData* currentData=(*it)->leafData().at(i).at(j);
		totalImportantNodes++;
		if (currentData!=NULL){
		    if(currentData->count()>=threshold&& currentData->revCount()>=threshold){
			//mNodesChecked[i][j]++;
			if (!currentData->isTrash()){
			    checkVariants(currentData);
			                        mNodesChecked[i][j]++;  
			  cout<<currentData->consensusFwd().length()<<endl;
			   for (int k=0; k<currentData->consensusFwd().length(); ++k){
				if (currentData->consensusFwd()[k]=='A'){
					mAs[i][j]++;
				}
				else if (currentData->consensusFwd()[k]=='C'){
					mCs[i][j]++;
				}
				else if (currentData->consensusFwd()[k]=='G'){
					mGs[i][j]++;
				}
				else if (currentData->consensusFwd()[k]=='T'){
					mTs[i][j]++;
				}
				else if (currentData->consensusFwd()[k]=='N'){
					mNs[i][j]++;
				}
			    }
			    vector <int> currentVariants = currentData->variants();
			    if (currentVariants.size()<10){
				while (!currentVariants.empty()){//************variants for each node consist of avector of ints, representing the hashes of each variant found trio*shift. of a 256 4*4*4*4 array made flat
				    int currentVariant = currentVariants.back();
				    currentVariants.pop_back();
				    mVariants[i][j][currentVariant]++;
				    mVariantsCount[i][j]++;

				}
				
			    }

			}

			else{cout<<"Trash"<<endl;}

		    }
		}
	    } 
	}
    }
    cout<<totalImportantNodes<<"=number of important nodes"<<endl;
}

void Trie::printTrieImportantOnly(Node* pCurrentNode, string barcode, int index){
//this function needs to be called before populate variants
    if ( pCurrentNode == NULL ){//if this is the first iteration, set current at root of trie
        pCurrentNode = mRootPointer;
        cout<<"Barcode Count"<<endl;
        barcode=string(mBarcodeLength, '\0');
    }
    else{//add the content of this node to the barcode
        barcode[index] = pCurrentNode->content();
        index++;
    }
    vector <Node*> children = pCurrentNode->children();
    if ( !children.empty() ){//go to next level of trie
        for (int i=0; i<children.size(); i++){
            pCurrentNode = children[i];
            printTrieImportantOnly(pCurrentNode, barcode, index);
        }
    }
    else if( !pCurrentNode->leafData().empty() ){//if we reach a leaf, print the count and variants
        ofstream summaryFile;
        summaryFile.open("/mnt/brick2/justin/SRR1613972_duplex/summaryIMPORTANT.txt", ios::app);
        for (int i=0; i<mNumberOfROIs; ++i){
            for (int j=0; j<mNumberOfPhases; ++j){
                LeafData* currentData= pCurrentNode->leafData()[i][j];
                if (currentData!=NULL && !currentData->isTrash() && mImportantNodes[i][j].find(pCurrentNode)!=mImportantNodes[i][j].end()){
                    summaryFile<<barcode<<" "<<mGenes[i]<<" phase "<<j<<endl;
                    summaryFile<<currentData->count()<<" reads"<<endl;
		    checkVariants(currentData); //this gets called twice... once in populate. there's pyobably a better way to do this
                    if (!currentData->variants().empty()){
                        for (int q=0; q<currentData->variants().size(); ++q){
                            summaryFile<<" "<<unhashVariants(currentData->variants()[q]).first<<" "<< unhashVariants(currentData->variants()[q]).second<<endl;
                        }
                    }
                }
            }
        summaryFile.close();
       }

       return;
    }
}


void Trie::printVariants(int threshold){
    cout<<"printing trie "<<endl;
    for (int i=0; i<mNumberOfROIs; ++i){
        for (int j=0; j<mNumberOfPhases; ++j){
            if (mVariantsCount[i][j]!=0){
                ostringstream os;
                os<<j;
                ostringstream os2;
                os2<<threshold;
                string filename= mGenes[i]+"_"+os.str()+"_thresh"+os2.str()+".txt";
                string matrixFilename = mGenes[i]+"_"+os.str()+"_thresh"+os2.str()+"_matrix.txt";
                ofstream outfile;
                ofstream matrixOutfile;
                outfile.open (filename.c_str());
                //matrixOutfile.open (matrixFilename.c_str());

                outfile<<"ROI: "<<mGenes[i]<<endl<<"Phase: "<<j<<endl<<"Total nodes checked: "<< mNodesChecked[i][j]<<endl<<"Total variants found: "<<mVariantsCount[i][j]<<endl;
                
		outfile<<mAs[i][j]<<" As"<<endl;
		outfile<<mCs[i][j]<<" Cs"<<endl;
		outfile<<mGs[i][j]<<" Gs"<<endl;
		outfile<<mTs[i][j]<<" Ts"<<endl;
		outfile<<mNs[i][j]<<" Ns"<<endl;
		/*
		map<int,int>::iterator it1;
                for (int l=0; l<5; ++l){//go through each base
                    for (int k = 0; k<mTargetLength[i]; ++k){
                        it1=mSubstitutions[i][j].find(k*5+l);
                        if (it1 == mSubstitutions[i][j].end()){
                            matrixOutfile<<left<<setw(15)<<setfill(' ')<<"0";   
                        }
                        else {
                            matrixOutfile<<left<<setw(15)<<setfill(' ')<<it1->second/(float)mSubstitutionsCount[i][j];
                        }
                    }
                    matrixOutfile<<endl;
                }*/
		for (int k=0; k<mVariants[i][j].size(); ++k){
		    if (mVariants[i][j][k]!=0){
           	         outfile<<unhashVariants(k).first<<" -> "<<unhashVariants(k).second<<" "<<mVariants[i][j][k]<<endl;
                    }
		}
                outfile.close();
            }
        }
    }  
}

void Trie::populateAndPrintVariants(){
    for (int i=0; i<mThresholdsOfImportance.size(); ++i){
	cout<<"populating"<<endl;
        populateVariants(mThresholdsOfImportance[i]);
        printVariants(mThresholdsOfImportance[i]);
    }
}
/* WHAT do we actually want for this? For a given ROI and a given PHase, the barcode count? OR the total count regardless?
int Trie::returnBarcodeCount(string barcode){
    Node* pCurrentNode = mRootPointer;
    int barcodeCount=0;
    for (int i = 0; i <= barcode.length(); i++){//go through the trie until you get to the end of the barcode
        Node* pNodeAtNextLevel = pCurrentNode->findChild(barcode[i]);
        if ( pNodeAtNextLevel == NULL ){
            if ( i == barcode.length() ){        
                barcodeCount = pCurrentNode->count();
            }
            return barcodeCount;
        }
        pCurrentNode = pNodeAtNextLevel;
    }
    return barcodeCount;
}
*/


void Trie::printTrie(Node* pCurrentNode, string barcode, int index){


    if ( pCurrentNode == NULL ){//if this is the first iteration, set current at root of trie
        pCurrentNode = mRootPointer;
        cout<<"Barcode Count"<<endl;
        barcode=string(mBarcodeLength, '\0');
    }
    else{//add the content of this node to the barcode
        barcode[index] = pCurrentNode->content();
        index++;
    }
    vector <Node*> children = pCurrentNode->children();
    if ( !children.empty() ){//go to next level of trie
        for (int i=0; i<children.size(); i++){
            pCurrentNode = children[i];
            printTrie(pCurrentNode, barcode, index);
        }
    }
    else if( !pCurrentNode->leafData().empty()){//if we reach a leaf, print the count and variants
        ofstream summaryFile;
        summaryFile.open("summary.txt", ios::app);
        summaryFile<<"barcode: "<<barcode<<endl;
        for (int i=0; i<mNumberOfROIs; ++i){
            summaryFile<<mGenes[i]<<endl;
            for (int j=0; j<mNumberOfPhases; ++j){
                LeafData* currentData= pCurrentNode->leafData()[i][j];
                if (currentData!=NULL){
                    summaryFile<<"phase "<<j<<endl;

                    if (!currentData->variants().empty()){
                        for (int q=0; q<currentData->variants().size(); ++q){
                            summaryFile<<" "<<unhashVariants(currentData->variants()[q]).first<<" "<< unhashVariants(currentData->variants()[q]).second<<endl;
                        }
                    }
                
                }
            }
        summaryFile.close();
       }

       return;
    }
}



/*and here. Maybe a barcode needs to have a total count at the node level
int Trie::returnMaxCount(int& max,Node* pCurrentNode ){
    if ( pCurrentNode==NULL ){ //if this is the first time the function is called, set current to root
        pCurrentNode=mRootPointer;
    }

    vector <Node*> children = pCurrentNode->children();

    if ( !children.empty() ){ //if there are still children, go one level down and recurse
        for (int i=0; i<children.size(); i++){
            pCurrentNode = children[i];
            returnMaxCount(max, pCurrentNode);
        }
    }
    else if ( pCurrentNode->count()>max ){//if the current count is greater than the max, reset the max
        max=pCurrentNode->count();
    }
    return max;
}

*/
