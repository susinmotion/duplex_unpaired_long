#ifndef VARIANT_H
#define VARIANT_H

#include <iostream>
#include <string>
using namespace std;

class Variant {
	int pos;
	char targ;
	char act;
	string trio;
	int hash;
	public:
	Variant(){};
	Variant(int p, char a, char t='\0', string tr="") {pos=p; targ=t; act=a; trio=tr; hash=hashVariants(p, t, a);}
	bool operator == (const Variant &ref) {
		if ( (this->pos == ref.getPos()) && (this->act==ref.getAct()) && (this->targ==ref.getTarg()) && (this->trio==ref.getTrio()) ){
			return true;
		}
		return false;
	}
	int hashVariants(int p, char t, char a){string nucleotides="ACGTN"; return p*20+nucleotides.substr(0,4).find(t)*5+nucleotides.find(a);}
	int getShiftHash(){return hash%20;}
	int getHash() const {return hash;}
	int getPos() const {return pos;}
	char getTarg()  const {return targ;}
	char getAct() const {return act;}
	string getTrio() const {return trio;}
};
#endif
