#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;
#include "TMath.h"

const char *FILE_TO_RECON = "analysis_results/ForNtrigStsRecon_.root";
const char *FILE_TO_RECON_MC = "150m.root";
bool DEBUG = false;
static const int    N_STN_MIN =   4;
static const int    N_STN_TOT =  19;
static const int    MATCH_HiSCORE = 1;
double dt_syst[N_STN_TOT];
TH1D *N_trig_sts = new TH1D("N_trig_sts","Number of stations triggered in event", N_STN_TOT - 1,  1,  N_STN_TOT);
TH1D *N_trig_stsMC = new TH1D("N_trig_stsMC","Number of stations triggered in event", N_STN_TOT - 1,  1,  N_STN_TOT);


int Draw(const char *FileName1 = FILE_TO_RECON, const char *FileName2 = FILE_TO_RECON_MC, int draw = 1)
{
    vector<Event> evts;
    ReadTreeOfEvts(FileName1, evts);
    for(int i = 0; i < (int) evts.size(); i++)
    {
        if(evts[i].NStations() >= N_STN_MIN && evts[i].Chi2() < 10){
            N_trig_sts -> Fill(evts[i].NStations());
        }
    }

    vector<Event> evtsMC;
    ReadTreeOfEvts(FileName2, evtsMC);
    for(int i = 0; i < (int) evtsMC.size(); i++)
    {
        if(evtsMC[i].NStations() >= N_STN_MIN){
            N_trig_stsMC -> Fill(evtsMC[i].NStations());
        }
    }

    if(draw)
    {
        auto c6 = new TCanvas("c6", "");
        N_trig_sts -> GetXaxis() -> SetTitle("number of stations");
        N_trig_sts -> Scale(1. / N_trig_sts -> GetEntries());
        N_trig_stsMC -> Scale(1. / N_trig_stsMC -> GetEntries());

        N_trig_sts -> SetLineWidth(2);

        //N_trig_sts -> Scale(1. / N_trig_sts -> GetEntries());
        TF1 *func = new TF1("func", "exp(-(x - 4) * 200. / [0]) * pow([1], x - 4) * [2]", 4, 10);
        func -> FixParameter(2, N_trig_sts -> GetMaximum());
        func -> SetParLimits(0, 1, 1000);
        func -> SetParLimits(1, 0, 1.);
        //N_trig_sts -> Fit("func","", "", 4, 10);

        //N_trig_stsMC -> Scale(N_trig_stsMC -> GetEntries());
        //N_trig_stsMC -> SetLineWidth(2);
        N_trig_stsMC -> SetLineColor(kRed);
        N_trig_sts -> SetMaximum(0.8);

        N_trig_sts -> Draw("P E SAME");
        N_trig_stsMC -> Draw("HIST SAME");
    }
    return 0;
}