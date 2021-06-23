#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;


const char *FILE_TO_DRAW = "analysis_results/Gr+Hisc221119.root";
int         N_STN_MIN = 2;
int         chanNumber = 4;
TH2D *points = new TH2D("points", "sky coordinates, x=ra, y=dec", 360, 0, 360, 90, -20, 70);
TH2D *points_up = new TH2D("points_up", "sky coordinates, x=ra, y=dec", 360, 0, 360, 90, -20, 70);
TH2D *points_up_down = new TH2D("points_up_down", "sky coordinates, x=ra, y=dec", 360, 0, 360, 90, -20, 70);
double pars[2][19];
//===============================================================================================//
int CalculateAverageDayPedestalPars(vector<EventAnalysis> &evts)
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
        TF1 *func = new TF1("func","gaus");// this variant is better to obtain RMS
        pedestal -> Fit(func,"Q");
        pars[0][stID] = func -> GetParameter(1);
        pars[1][stID] = func -> GetParameter(2);
        delete pedestal;
    }
    return 0;
}
//===============================================================================================//
int SetNewPedestalParsChan4(vector<EventAnalysis> &evts)
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
//===============================================================================================//
int DrawSkyPoints(const char *fileName = FILE_TO_DRAW, int draw = 1)
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
    CalculateAverageDayPedestalPars(evts);
    //for(int i=0; i < 19; i++) cout << pars[1][i]<< endl;
    SetNewPedestalParsChan4(evts);
   // cout << evts[1].MuStn(0).ID() << " " << evts[1].EvtPedestalRMS(evts[1].MuStn(0).ID()-1)[3] << endl;

    for(int i = 0; i < (int)evts.size(); i++) evts[i].FindSignal();

    for(int i = 0; i < (int)evts.size(); i++) if(evts[i].MatchHiSC() == 1 && evts[i].NStations() >= N_STN_MIN) {
        points->Fill(evts[i].RaHiSC(), evts[i].DecHiSC());
        if(evts[i].SigStart(evts[i].MuStn(0).ID()-1)[chanNumber] != -1000) points_up_down->Fill(evts[i].RaHiSC(), evts[i].DecHiSC());
        if(evts[i].SigStart(evts[i].MuStn(0).ID()-1)[chanNumber] == -1000) points_up->Fill(evts[i].RaHiSC(), evts[i].DecHiSC());
        }
    if(draw)
    {
        auto c4 = new TCanvas("c4", "all");
        points->Draw("CONT4Z");
        auto c5 = new TCanvas("c5", "up");
        points_up->Draw("CONT4Z");
        auto c6 = new TCanvas("c6", "up-down");
        points_up_down->Draw("CONT4Z");
    }
    return 0;
}
//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
int main1()
{
    DrawSkyPoints("analysis_results/Gr+Hisc221119_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc241119_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc261119_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc271119_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc291119_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc241219_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc251219_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc291219_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc301219_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc020120_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc130220_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc150220_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc160220_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc200220_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc210220_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc250220_newFolder.root",0);
    DrawSkyPoints("analysis_results/Gr+Hisc270220_newFolder.root");
    return 0;
}