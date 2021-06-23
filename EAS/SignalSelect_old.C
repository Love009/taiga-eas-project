#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;

const int TIME_WINDOW = 5; // *5 ns
vector<double> MeanCH0;
vector<double> MeanCH2;
vector<double> SigmaCH0;
vector<double> SigmaCH2;
TH1D *signalStart = new TH1D("signalStart","signalStart", 200,  0,  200); //hist for linear fit
const char *INITIALIZATION_FILE_NAME = "analysis_results/pedestals/pedestal291219.txt";
//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
//Initualize station channel pedestals
void Initialize()
{
    double meanCH0, meanCH2, sigmaCH0, sigmaCH2;
    string  title, rest;
    int i = 1;
    ifstream in(INITIALIZATION_FILE_NAME);
    getline(in, title);
    while (1) {
        in >> rest >> meanCH0 >> sigmaCH0 >> meanCH2 >> sigmaCH2;
        if (!in.good()) break;
        MeanCH0.push_back(meanCH0);
        MeanCH2.push_back(meanCH2);
        SigmaCH0.push_back(sigmaCH0);
        SigmaCH2.push_back(sigmaCH2);
        i++;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Find a signal start using linear fit (the least square method)
int FindSignal(vector<int> ADCsig, int stID, int chanID, double sigma, int threshold, int lowLimit)
{
    double sigTime = -1., flag = 0;
    for(int k = 0; k < (int) ADCsig.size(); k++)
    {
        double mean = MeanCH0[stID - 1];
        if(chanID == 2) mean = MeanCH2[stID - 1];
        int delta = ADCsig[k] -  mean;
        if(delta <= lowLimit * sigma) 
        {
            cout << "Signal ampl is lower than " << lowLimit * sigma << endl;
            break;
        }
        else if(flag >= TIME_WINDOW)
        {
            sigTime = k - 1 - TIME_WINDOW;
            break;
        }
        else if(delta >= threshold * sigma) flag++;
        else flag = 0;
    }
    return (sigTime * 5);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CorrectedTree(const char * FileName = "analysis_results/Gr+Hisc291219.root", const char * FileToWrite = "analysis_results/Gr+Hisc291219_Generation2.root", int threshold = 5, int lowLimit = -10)
{
    vector<Event> evts;
    ReadTreeOfEvts(FileName, evts);
    cout << evts.size() << " events before" << endl;
    Initialize();
    int nSigTot = 0, nSig = 0, sigTime0/*ns*/, sigTime2/*ns*/, counter = 0, counter_init = 0;
    
    for(int i = 0; i < (int) evts.size(); i++)
    {
        for(int j = 0; j < evts[i].NStations(); j++)
        {
            int stID = evts[i].MuStn(j).ID();
            counter_init++;
            sigTime0 = FindSignal(evts[i].MuStn(j).ADC0(), stID, 0, SigmaCH0[stID - 1], threshold, lowLimit);
            sigTime2 = FindSignal(evts[i].MuStn(j).ADC2(), stID, 2, SigmaCH2[stID - 1], threshold, lowLimit);
            if(sigTime0 > 0 && sigTime2 > 0)  
            {
                int sigTime = sigTime0;
                if(sigTime2 < sigTime) sigTime = sigTime2;
                double newT = evts[i].MuStn(j).T() + (double)sigTime;
                evts[i].MuStn(j).SetT(newT);
                counter++;
            }
            else {evts.erase(evts.begin() + i); i--;}
        }
    }
    
    cout << "Selected " << counter << " of " << counter_init << " signals" << endl;
    cout << evts.size() << " events after"<< endl;
//Write new root-file   
    TFile f(FileToWrite,"recreate");
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
