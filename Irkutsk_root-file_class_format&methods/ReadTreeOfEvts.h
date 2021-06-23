#include <stdlib.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include "TVector.h"
using namespace std;
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "AutoDict_Event_cxx.so"
#else
#include "TunkaIrkutsk.h"
#endif

//Read vector of Events from a root-file and put it in evts
void ReadTreeOfEvts(const char *FileName, vector<Event> &evts)
{
    int stsevts = 0;
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
        if(event->NStations() == 4) stsevts++;
        event->Reset();
    }
    cout << nEvTot << " events" << endl;
    cout << stsevts << " 4 st events" << endl;
    f->Close();
}


