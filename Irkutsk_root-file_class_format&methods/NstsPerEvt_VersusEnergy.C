//Number of stations dependence on EAS energy

#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;
const char *FILE_NAME = "analysis_results/Gr+Hisc130220_newFolder.root";
TH2D * N_E = new TH2D("N_E", "N(E), hist", 19, 1, 20, 2048*7, 0, 2048*7);
TH1D * E = new TH1D("E", "E, hist", 2048*7, 0, 2048*7);
TProfile * pr_N_E = new TProfile("pr_N_E","profile N_E",19, 1, 20, 0, 2048*7);
double pars[2][19];
//===============================================================================================//
int CalculateAverageDayPedestalPars(vector<EventAnalysis> &evts, int chanNumber)
{
    for(int stID=0; stID < 19; stID++)
    {
        TH1D *pedestal = new TH1D("pedestal", "Pedestal", 80, 2030, 2070);
        for (int k = 0; k < (int) evts.size(); k++)
        {
            EventAnalysis evt = evts[k];
            for(int l = 0; l < evt.NStations(); l++)
            {
                int stID_cur = evt.MuStn(l).ID() - 1;
                if(stID_cur == stID) pedestal->Fill(evt.EvtPedestalMean(stID_cur)[chanNumber]);
            }
        }
        TF1 *func = new TF1("func","gaus");// this variant is better
        pedestal -> Fit(func,"Q");
        pars[0][stID] = func -> GetParameter(1);
        pars[1][stID] = func -> GetParameter(2);
        delete pedestal;
    }
    return 0;
}
//===============================================================================================//
int SetNewPedestalParsChan4(vector<EventAnalysis> &evts, int chanNumber)
{
    for (int k = 0; k < (int) evts.size(); k++)
    {
        evts[k].SetEvtPedestalMean(pars[0], chanNumber);
        evts[k].SetEvtPedestalRMS(pars[1], chanNumber);
    }
    return 0;
}
//===============================================================================================//
//===============================================================================================//
int NstsPerEvt_VersusEnergy(const char *fileName = FILE_NAME, int draw = 1)//ground stations are considered
{
    vector<Event> evts0;
    ReadTreeOfEvts(fileName, evts0);
    vector<EventAnalysis> evts;
    for(int i = 0; i < (int)evts0.size(); i++)
    {
        EventAnalysis evt(evts0[i]);
        evt.CalculatePedestalPars();
        evts.push_back(evt);
    }

    for(int chanNumber = 0; chanNumber < 4; chanNumber++)
    {
        CalculateAverageDayPedestalPars(evts, chanNumber);
        SetNewPedestalParsChan4(evts, chanNumber);
    }
    double mean = 0.;
    int norm;
    for(int i = 0; i < (int)evts.size(); i++)
    {
        mean = 0.;
        evts[i].FindSignal();
        evts[i].FindSigMax();
        norm = evts[i].NStations();

        if(evts[i].NStations() > 1 && evts[i].MatchHiSC() == 1)
        {
            if(abs(evts[i].XYHiSC()[0]) <= 400 && abs(evts[i].XYHiSC()[1]) <= 400) {


                for (int j = 0; j < evts[i].NStations(); j++) {
                    int stationID = evts[i].MuStn(j).ID() - 1;
                    if (evts[i].SigMax(stationID)[0] != -1000 && evts[i].SigMax(stationID)[2] != -1000) {
                        mean += (evts[i].SigMax(stationID)[0] + evts[i].SigMax(stationID)[2]) / 2.;
                    } else if (evts[i].SigMax(stationID)[1] != -1000 &&
                               evts[i].SigMax(evts[i].MuStn(j).ID())[3] != -1000) {
                        mean += (evts[i].SigMax(stationID)[1] + evts[i].SigMax(stationID)[3]) * 10 / 2.;
                    } else norm--;
                }
                // in ideal case norm = evts[i].NStations()
                if (mean != 0 && norm > 1) {
                    N_E->Fill(norm, mean);
                    E->Fill(mean);
                    pr_N_E->Fill(norm, mean);
                }
            }}
    }
    if(draw)
    {
        auto c4 = new TCanvas("c4", "");
        E -> Draw();
        auto c5 = new TCanvas("c5", "");
        N_E->GetXaxis()->SetTitle("N stations");
        N_E->GetYaxis()->SetTitle("E");
        // N_E -> SetMarkerStyle(6);
       // N_E -> SetMarkerSize(1);
        N_E -> Draw("BOX");
        auto c6 = new TCanvas("c6", "");
        pr_N_E -> Draw("HIST TEXT0");
        pr_N_E -> Draw("B SAME");
    }
    return 0;
}
int main1()
{
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc221119_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc241119_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc261119_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc271119_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc291119_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc201219_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc241219_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc251219_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc291219_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc301219_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc020120_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc130220_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc150220_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc160220_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc200220_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc210220_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc250220_newFolder.root", 0);
    NstsPerEvt_VersusEnergy("analysis_results/Gr+Hisc270220_newFolder.root");
    return 0;
}
