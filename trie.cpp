#include "trie.h"
#include "node.h"
#include "variants.h"
#include "variant.h"
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
#include <algorithm>
using namespace std;

Node* Trie::pRootPointer(){
    return mRootPointer;
}

void Trie::addBarcode(int ROINumber, int phase, string barcode, string sequence, string target, string rev){
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
	    if (rev=="rev"){
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
    mSuperVariantsCount=vector <vector<int> >(mNumberOfROIs, vector<int>(mNumberOfPhases,  0));
    mNodesChecked= vector <vector<int> >(mNumberOfROIs, vector<int>(mNumberOfPhases,  0));
    mVariants=vector <vector <vector<int> > >(mNumberOfROIs, vector<vector< int> >(mNumberOfPhases, vector<int>(472500, 0)));
    Variant* emptyVariant=new Variant();
    mSuperVariants=vector <vector <vector<int> > >(mNumberOfROIs, vector<vector<int> >(mNumberOfPhases, vector<int>(427500, 0)));           
    mSuperShifts=vector<vector <vector<int> > >(mNumberOfROIs, vector<vector <int> >(mNumberOfPhases, vector<int>(20,0)) );
    mShifts=vector<vector <vector<int> > >(mNumberOfROIs, vector<vector <int> >(mNumberOfPhases, vector<int>(20,0)) );
    cout<<mNumberOfROIs<< " rois"<<endl;
    cout<<mNumberOfPhases<<" phases"<<endl;
    int totalImportantNodes=0;
    for (int i=0; i<mNumberOfROIs; ++i){
	for (int j=0; j<mNumberOfPhases; ++j){
	    set <Node*> currentImportantNodes=mImportantNodes[i][j];
	    cout<<currentImportantNodes.size()<<" important nodes"<<endl;
            ofstream dupes;
            dupes.open("dupes.txt");
	    for (set <Node*>::iterator it=currentImportantNodes.begin(); it!=currentImportantNodes.end(); ++it){
		LeafData* currentData=(*it)->leafData().at(i).at(j);
		
		totalImportantNodes++;
		if ( (totalImportantNodes%10000)==0){
			cout<<totalImportantNodes<<endl;
		}
		if (currentData!=NULL){
	              if( (currentData->count()-currentData->revCount()>=threshold)&& (currentData->revCount()>=threshold) ){
			mNodesChecked[i][j]++;
			if (!currentData->isTrash()){
			    vector <Variant*> currentFwdVariants=bowtieCheckVariants(currentData->consensusFwd(), "mus", currentData);
			    vector <Variant*> currentRevVariants=bowtieCheckVariants(currentData->consensusRev(), "mus");
			    if (currentFwdVariants.size()!=0 or currentRevVariants.size()!=0){
				sort(currentFwdVariants.begin(), currentFwdVariants.end());
				sort(currentRevVariants.begin(), currentRevVariants.end());
				vector <Variant*> currentSuperVariants;
				vector<int>hashes;
				for (vector<Variant*>::iterator it=currentFwdVariants.begin(); it!=currentFwdVariants.end();++it){
					vector<Variant*>::iterator it1=currentRevVariants.begin();
					for (; it1!=currentRevVariants.end(); ++it1){
						if ( (*it)->getPos()==(*it1)->getPos()){
							if( (*it)->getHash()!=(*it1)->getHash()){
								cout<<"dupe!"<<endl;
								dupes<<(*it)->getPos()<<" "<<(*it)->getTarg()<<" -> "<<(*it)->getAct()<<" : "<<(*it1)->getPos()<<" "<<(*it1)->getTarg()<<" -> "<<(*it1)->getAct()<<endl;
							}
							else{   

								mSuperVariants[i][j][(*it)->getHash()]++;
		                                                mSuperShifts[i][j][(*it)->getShiftHash()]++;
                		                                mSuperVariantsCount[i][j]++;
							}
							currentRevVariants.erase(it1);
							break;
						}
					}
					if (it1==currentRevVariants.end()){
						mVariants[i][j][(*it)->getHash()]++;
						mShifts[i][j][(*it)->getShiftHash()]++;
						mVariantsCount[i][j]++;
					}
				}
				for (vector<Variant*>::iterator it=currentRevVariants.begin(); it!=currentRevVariants.end(); ++it){
					mVariants[i][j][(*it)->getHash()]++;
                                        mShifts[i][j][(*it)->getShiftHash()]++;
                                        mVariantsCount[i][j]++;
				}
/*
			    	set_intersection(currentFwdVariants.begin(), currentFwdVariants.end(), currentRevVariants.begin(), currentRevVariants.end(), back_inserter(currentSuperVariants));
			    	vector <Variant*> currentStrandDifs;
			    	set_symmetric_difference(currentFwdVariants.begin(), currentFwdVariants.end(), currentRevVariants.begin(), currentRevVariants.end(), back_inserter(currentStrandDifs));
			    	if (currentSuperVariants.size()<10){
					while (!currentSuperVariants.empty()){
						Variant* currentSuperVariant=currentSuperVariants.back();
						currentSuperVariants.pop_back();
						mSuperVariants[i][j][currentSuperVariant->getHash()]++;
						mSuperShifts[i][j][currentSuperVariant->getShiftHash()]++;
						mSuperVariantsCount[i][j]++;
					}
				}		
				    if (currentStrandDifs.size()<10){
					ofstream dupes;
					dupes.open("dupes.txt", ofstream::app);
					for (int k=0; k<currentStrandDifs.size(); ++k){
						Variant* currentDif=currentStrandDifs.at(k);
						int l=k+1;	
					if ( (l<currentStrandDifs.size()) && (currentDif->getPos()==currentStrandDifs.at(l)->getPos()) ){
							cout<<"dupes!"<<endl;
							k++;
							dupes<<currentDif->getPos()<<" "<<currentDif->getTarg()<<" -> "<<currentDif->getAct()<< " : "<<currentStrandDifs.at(l)->getPos()<<" "<<currentStrandDifs.at(l)->getTarg()<<" -> "<<currentStrandDifs.at(l)->getAct()<<endl;
						}
						else{
						mVariants[i][j][currentDif->getHash()]++;
						mShifts[i][j][currentDif->getShiftHash()]++;
				   		mVariantsCount[i][j]++;
						}
					}
					dupes.close();

			    	    }
				}
*/	
				}
			}	
		else{cout<<"Trash"<<endl;}

		    }
		}dupes.close(); 
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
        summaryFile.open("summaryIMPORTANT.txt", ios::app);
        for (int i=0; i<mNumberOfROIs; ++i){
            for (int j=0; j<mNumberOfPhases; ++j){
                LeafData* currentData= pCurrentNode->leafData()[i][j];
                if (currentData!=NULL && !currentData->isTrash() && mImportantNodes[i][j].find(pCurrentNode)!=mImportantNodes[i][j].end()){
                    summaryFile<<barcode<<" "<<mGenes[i]<<" phase "<<j<<endl;
                    summaryFile<<currentData->count()<<" reads"<<endl;
                    if (!currentData->variants().empty()){
                        for (int q=0; q<currentData->variants().size(); ++q){
                            //summaryFile<<" "<<unhashVariants(currentData->variants()[q]).first<<" "<< unhashVariants(currentData->variants()[q]).second<<endl;
                        	summaryFile<<" "<<currentData->variants()[q]<<endl;//this is counting on my overloading << and also making this a vector of variant pointers not a vector of ints
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
            if (mVariantsCount[i][j]!=0 || mSuperVariantsCount[i][j]!=0){
                ostringstream os;
                os<<j;
                ostringstream os2;
                os2<<threshold;
                string filename= mGenes[i]+"_"+os.str()+"_thresh"+os2.str()+".txt";
                string matrixFilename = mGenes[i]+"_"+os.str()+"_thresh"+os2.str()+"_matrix.txt";
		string supermatrixFilename=mGenes[i]+"_"+os.str()+"_thresh"+os2.str()+"_supermatrix.txt";
                string superFilename=mGenes[i]+"_"+os.str()+"_thresh"+os2.str()+"_super.txt";
                ofstream outfile;
                ofstream matrixOutfile;
		ofstream superOutfile;
                ofstream supermatrixOutfile;
		outfile.open (filename.c_str());
                matrixOutfile.open (matrixFilename.c_str());
		supermatrixOutfile.open (supermatrixFilename.c_str());
		superOutfile.open(superFilename.c_str());
                outfile<<"ROI: "<<mGenes[i]<<endl<<"Phase: "<<j<<endl<<"Total nodes checked: "<< mNodesChecked[i][j]<<endl<<"Total variants found: "<<mVariantsCount[i][j]<<endl;
                superOutfile<<"ROI: "<<mGenes[i]<<endl<<"Phase: "<<j<<endl<<"Total nodes checked: "<< mNodesChecked[i][j]<<endl<<"Total variants found: "<<mSuperVariantsCount[i][j]<<endl;

		map<int,int>::iterator it1;
		
		for (int l=0; l<4; ++l){
                        for (int k=0; k<5; ++k){
				supermatrixOutfile<<left<<setw(15)<<setfill(' ')<<mSuperShifts[i][j][l*5+k]/float(mSuperVariantsCount[i][j]);
                                matrixOutfile<<left<<setw(15)<<setfill(' ')<<mShifts[i][j][l*5+k]/float(mVariantsCount[i][j]);
                        }
                       supermatrixOutfile<<endl;
			 matrixOutfile<<endl;
                }
               for (int l=0; l<4; ++l){
                        for (int k=0; k<5; ++k){
                                supermatrixOutfile<<"ACGTN"[l]<<" -> "<<"ACGTN"[k]<<" "<<mSuperShifts[i][j][l*5+k]<<endl;
                        	matrixOutfile<<"ACGTN"[l]<<" -> "<<"ACGTN"[k]<<" "<<mShifts[i][j][l*5+k]<<endl;
			}
                }
                matrixOutfile.close();
		supermatrixOutfile.close();

		/*
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
               			string nucleotides="ACGTN";
                                outfile<<k/20<<" "<<nucleotides[(k%20)/5]<<" -> "<<nucleotides.substr(k%5, 1)<<" "<<mVariants[i][j][k]<<endl;    
		    }
		}
                outfile.close();
		for (int k=0; k<mSuperVariants[i][j].size();++k){
			if (mSuperVariants[i][j][k]!=0){
				string nucleotides="ACGTN";
				superOutfile<<k/20<<" "<<nucleotides[(k%20)/5]<<" -> "<<nucleotides.substr(k%5, 1)<<" "<<mSuperVariants[i][j][k]<<endl;
			}
		}
		
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
                        //summaryFile<<" "<<unhashVariants(currentData->variants()[q]).first<<" "<< unhashVariants(currentData->variants()[q]).second<<endl;
                                summaryFile<<" "<<currentData->variants()[q]<<endl;//this is counting on my overloading << and also making this a vector of variant pointers not a vector of ints
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
