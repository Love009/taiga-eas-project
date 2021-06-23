#include "common.h"
using namespace std;
#include <stdlib.h>
#include <bitset>
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "AutoDict_Event_cxx.so"
#else
#include "Tunka.h"
#endif

TH1D *pedestal0 = new TH1D("pedestal0","Pedestal channel 0", 80, 2030, 2070);
TH1D *pedestal2 = new TH1D("pedestal2","Pedestal channel 2", 80, 2030, 2070);
double pars0[2];//mean and sigma of CH0
double pars2[2];//mean and sigma of CH2
const char *INPUT_FILE_NAME = "analysis_results/Gr271119.root";
const char *OUTPUT_FILE_NAME = "analysis_results/pedestals/pedestal271119.txt";
static const int    N_STN_TOT =  19;
//=========================================================================================================
//Find pedestal for st with stID, for 0 and 2 channels and draw two histograms
int Pedestal(const char * FileName = INPUT_FILE_NAME, int stID = 1)
{
    int counter = 0;
    vector<Event> evts;
    TFile *f = new TFile(FileName);
    TTree *t4 = (TTree*)f->Get("Muon");
    Event *event = new Event();
    TBranch *br = t4->GetBranch("events");
    br ->SetAddress(&event);
    t4 -> SetAutoSave(1000000);
    long int nEvTot = t4 -> GetEntries();
    for(long int i = 0; i < nEvTot; i++)
    {
       br ->GetEntry(i);
       evts.push_back(*event);
       event->Reset();
    }
    f->Close();
    
    pedestal0->Reset();
    pedestal2->Reset();
    for(int i = 0; i < (int) evts.size(); i++)
    {
        for(int j = 0; j < evts[i].NStations(); j++)
        {
            if(evts[i].MuStn(j).ID() == stID)
            {
                double mean0 = 0, mean2 = 0;
                int norm0 = 20, norm2 = 20;//for calc of the average
                for(int k = 0; k < 20; k++)
                {
                    if(evts[i].MuStn(j).ADC0()[k] < 2048 + 50 && evts[i].MuStn(j).ADC0()[k] > 2048 - 50) mean0 += evts[i].MuStn(j).ADC0()[k];
                    else norm0--;
                    if(evts[i].MuStn(j).ADC2()[k] < 2048 + 50 && evts[i].MuStn(j).ADC2()[k] > 2048 - 50)mean2 += evts[i].MuStn(j).ADC2()[k];
                    else norm2--;
                    counter++;
                }
                pedestal0 -> Fill(mean0 / norm0);
                pedestal2 -> Fill(mean2 / norm2);
            }
        }
    }    
    TF1 *func = new TF1("func","gaus");
    pedestal0 -> Fit(func,"Q");
    pars0[0] = func -> GetParameter(1);
    pars0[1] = func -> GetParameter(2);   
    pedestal2 -> Fit(func,"Q");
    pars2[0] = func -> GetParameter(1);
    pars2[1] = func -> GetParameter(2);    
    //Drawing
    /*TCanvas *c1 = new TCanvas("c1", "Pedestal channel 0");
    pedestal0 -> Draw();
    TCanvas *c2 = new TCanvas("c2", "Pedestal channel 2");
    pedestal2 -> Draw();*/
    return 0;
}

int Pedestals(const char * FileName = INPUT_FILE_NAME, const char * FileToWrite = OUTPUT_FILE_NAME)
{
    ofstream out(FileToWrite);
    out << "stID\tmeanCH0\tsigmaCH0\tmeanCH2\tsigmaCH2" << "\n";
    for(int i = 1; i <= N_STN_TOT; i++)
    {
        Pedestal(FileName, i);
        out << i << std::setprecision(6) << "\t" << pars0[0] << "\t"  << std::setprecision(3) << pars0[1] << "\t" << std::setprecision(6) << pars2[0] << "\t" << std::setprecision(3) << pars2[1] << "\t" << endl;        
    }
    out.close();
    return 0;    
}
