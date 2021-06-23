#include "ReadTreeOfEvts.h"

//Drawing of oscillogramm
int DrawOsc_0_2(const char * FileName = "analysis_results/Gr+Hisc261119_Generation3.root", int stID = 2)
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
    TCanvas *c1 = new TCanvas("c1", "Signal from channel");
    TCanvas *c2 = new TCanvas("c2", "Signal from channel ");
    for(int i = 0; i < 10000 /* (int) evts.size()*/; i++)
    {
        TH1I *sig = new TH1I("sig","Signal from channel 0", evts[i].MuStn(0).ADC0().size(), 0, evts[i].MuStn(0).ADC0().size());
        TH1I *sig2 = new TH1I("sig2","Signal from channel 2", evts[i].MuStn(0).ADC2().size(), 0, evts[i].MuStn(0).ADC2().size());
        sig ->GetYaxis()->SetRangeUser(2000, 5000);
        sig2 ->GetYaxis()->SetRangeUser(2000, 5000);
        for(int j = 0; j < evts[i].NStations(); j++)
        {
            if(evts[i].MuStn(j).ID() == stID)
            {
                for(int k = 0; k < (int) evts[i].MuStn(j).ADC0().size(); k++)
                {
                    sig -> SetBinContent((k + 1), evts[i].MuStn(j).ADC0()[k]);
                    sig2 -> SetBinContent((k + 1), evts[i].MuStn(j).ADC2()[k]);
                    counter++;
                }
                c1->cd(1);
                if(i != 0)    sig->Draw("SAME");
                else sig -> Draw();
                c2 -> cd(1);
                if(i != 0)    sig2->Draw("SAME");
                else sig2 -> Draw();
            }
        }
    }
    cout << counter << endl;
    int counter2 = 0;
    for(int i = 0; i < (int) evts.size(); i++)
    {
      for(int j = 0; j < evts[i].NStations(); j++)  
      {
        if(evts[i].MuStn(j).ID() == 1) counter2++;  
      }
    }
    cout << counter2 * 200 << endl;
    return 0;
}
//Drawing of oscillogramm
int DrawOsc_1_3(const char * FileName = "analysis_results/Gr291119.01.root", int stID = 2)
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
    //TH1I *sig = new TH1I("sig","Signal from channel 1", evts[0].MuStn(0).ADC1().size(), 0, evts[0].MuStn(0).ADC1().size());
    TCanvas *c1 = new TCanvas("c1", "Signal from channel");
    TCanvas *c2 = new TCanvas("c2", "Signal from channel ");
    for(int i = 0; i < 10000 /* (int) evts.size()*/; i++)
    {
        //delete sig;
        TH1I *sig = new TH1I("sig","Signal from channel 1", evts[i].MuStn(0).ADC1().size(), 0, evts[i].MuStn(0).ADC1().size());
        TH1I *sig2 = new TH1I("sig2","Signal from channel 3", evts[i].MuStn(0).ADC3().size(), 0, evts[i].MuStn(0).ADC3().size());
        sig ->GetYaxis()->SetRangeUser(2000, 10000);
        sig2 ->GetYaxis()->SetRangeUser(2000, 10000);
        for(int j = 0; j < evts[i].NStations(); j++)
        {
            if(evts[i].MuStn(j).ID() == stID)
            {
                for(int k = 0; k < (int) evts[i].MuStn(j).ADC1().size(); k++)
                {
                    sig -> SetBinContent((k + 1), evts[i].MuStn(j).ADC1()[k]);
                    sig2 -> SetBinContent((k + 1), evts[i].MuStn(j).ADC3()[k]);
                    counter++;
                }
                c1->cd(1);
                if(i != 0)    sig->Draw("SAME");
                else sig -> Draw();
                c2 -> cd(1);
                if(i != 0)    sig2->Draw("SAME");
                else sig2 -> Draw();
            }
        }
    }
    cout << counter << endl;
    int counter2 = 0;
    for(int i = 0; i < (int) evts.size(); i++)
    {
        for(int j = 0; j < evts[i].NStations(); j++)
        {
            if(evts[i].MuStn(j).ID() == 1) counter2++;
        }
    }
    cout << counter2 * 200 << endl;
    return 0;
}
