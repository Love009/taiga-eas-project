#include "common.h"
using namespace std;
#include <stdlib.h>
#include <bitset>
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "AutoDict_Event_cxx.so"
#else
#include "Tunka.h"
#endif

//Read vector of Events from a root-file and put it in evts
void ReadTreeOfEvts(const char *FileName, vector<Event> &evts)
{
    TFile *f = new TFile(FileName);
    TTree *t4 = (TTree*)f->Get("Muon");
    Event *event = new Event();
    TBranch *br = t4->GetBranch("events");
    br ->SetAddress(&event);
    long int nEvTot = t4 -> GetEntries();
    for(long int i = 0; i < nEvTot; i++)
    {
       br ->GetEntry(i);
       evts.push_back(*event);
       event->Reset();
    }
    cout << nEvTot << " events" << endl;
    f->Close();
}
