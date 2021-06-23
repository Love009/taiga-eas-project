#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;
#include "TMath.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int WriteToFile(int ID, vector<Event> &evts, const char* FileName)
{
    string fileToWrite(FileName);
    fileToWrite = fileToWrite.substr(0, (int)strlen(FileName) - 5);
    fileToWrite += "_"+to_string(ID)+".root";
    TFile f(fileToWrite.c_str(),"recreate");
    TTree Muon("Muon","");
    gROOT->ProcessLine("#include <vector>");
    Event *event = new Event();
    Muon.Branch("events", &event);
    Muon.SetAutoSave(1000000);
    for(int i = 0; i < (int)evts.size(); i++)
    {
        event = & evts[i];
        Muon.Fill();
    }
    f.Write();
    f.Close();
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Split(const char *FileName)
{
    vector<Event> evts, evts_[3];
    ReadTreeOfEvts(FileName, evts);
    int counter2 = 0;
    for(int i=0; i < (int) evts.size(); i++)
    {
        if(evts[i].NStations() == 2 && evts[i].MatchHiSC() == 1) counter2++;
    }
    cout << counter2 << endl;
    for(int i = 0; i < 8000; i++)
    {
       evts_[0].push_back(evts[i]);
    }
    WriteToFile(1,evts_[0],FileName);
    for(int i = 8000; i < 16000; i++)
    {
        evts_[1].push_back(evts[i]);
    }
    WriteToFile(2,evts_[1],FileName);
    for(int i = 16000; i < (int) evts.size(); i++)
    {
        evts_[2].push_back(evts[i]);
    }
    WriteToFile(3,evts_[2],FileName);
    return 0;
}