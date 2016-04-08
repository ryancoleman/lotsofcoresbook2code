/* Code generator for QPhiX Library */

#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <typeinfo>
#include <string>
#include <map>

using namespace std;

#include "dslash.h"

void generate_code(void);

// Merge L2 prefetches with another instruction stream
void mergeIvectorWithL2Prefetches(InstVector& ivector, InstVector& l2prefs)
{
    if( l2prefs.size() == 0 ) {
        cout << "No L2 Prefetches. Returning ivector unchanged" << endl;
    }
    else {

        if( ivector.size() == 0 ) {
            cout << "No actual instructions to merge. Returning " << endl;
            return;
        }


        int ivector_size = ivector.size();
        vector<Instruction*>::iterator it;

        // Skip past first declarations
        it=ivector.begin();

        while ( (*it)->numDeclarations() > 0 ) {
            it++;
            ivector_size--;
        }

        int n_prefs = l2prefs.size();

        //cout << "After declarations Ivector size is " << ivector_size << endl;
        //cout << "PrefetchL2 size is " << n_prefs << endl;

        int spacing_factor = ivector_size / n_prefs;
        //cout << "Spacing Factor is  " << spacing_factor << endl;
        int pref=0;

        for(int i=0; i < n_prefs; i++) {
            it = ivector.insert(it, l2prefs[pref]);
            pref++;
            it++;
            int j=spacing_factor;

            while (j > 0) {
                it++;
                j--;
            }
        }
    }

    std::map<string,string> offslist;

    for(int i=0; i < ivector.size(); i++) {
        Instruction *inst = ivector[i];

        if ( inst->hasAddress() ) {
            MemRefInstruction* mr = dynamic_cast< MemRefInstruction* >(inst);

            if(mr->hasGSAddress()) {
                const GatherAddress* ga = dynamic_cast<const GatherAddress *>(mr->getAddress());
                string offs = ga->getOffsets(false);
                string voffs = ga->getOffsets ();

                if(offslist.find(offs) == offslist.end()) {
                    offslist[offs] = voffs;
                }
            }
        }
    }

    for ( map<string,string>::iterator  it = offslist.begin(); it != offslist.end(); ++it ) {
        ivector.insert(ivector.begin(), new DeclareOffsets(it->first, it->second ));
    }
}


// Dump an instruction stream into a file
void dumpIVector(InstVector& ivector, string filename)
{
    ofstream outfile(filename.c_str());

    for(int i=0; i < ivector.size(); i++) {
        outfile << ivector[i]->serialize() << endl;
    }

    outfile.close();
}

int main(int argc, char *argv[])
{
    generate_code();
    return 0;
}
