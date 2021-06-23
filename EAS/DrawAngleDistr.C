#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;
#include "TMath.h"

const char *FILE_TO_RECON = "analysis_results/ForNtrigStsRecon_.root";
const char *FILE_TO_RECON_MC = "150m_angle.root";
bool DEBUG = false;
static const int    N_STN_MIN =   3;
static const int    N_STN_TOT =  19;
static const int    MATCH_HiSCORE = 1;
double dt_syst[N_STN_TOT];
TH1D *thetaDistrMC = new TH1D("thetaDistrMC", "theta, degrees", 190,  -5,  180);
TH1D *phiDistrMC = new TH1D("phiDistrMC", "phi, degrees", 370 / 4.,  -5,  365);
TH1D *thetaDistrHiSC = new TH1D("thetaDistrHiSC", "theta, degrees", 190,  -5,  180);
TH1D *phiDistrHiSC = new TH1D("phiDistrHiSC", "phi, degrees", 370 / 4.,  -5,  365);
TH2D *points = new TH2D("points", "phi sinHiscTheta, theta", 360, -180, 180, 90, -90, 90);

int Draw(const char *FileName1 = FILE_TO_RECON, const char *FileName2 = FILE_TO_RECON_MC, int draw = 1)
{
    vector<Event> evts;
    ReadTreeOfEvts(FileName1, evts);
    double HiSC_the = 25. / (180. / TMath::Pi());
    for(int i = 0; i < (int) evts.size(); i++)
    {
        if(evts[i].NStations() >= N_STN_MIN && evts[i].MatchHiSC() == 1){
            thetaDistrHiSC->Fill(180. / TMath::Pi() * evts[i].ThetaGen());
            phiDistrHiSC->Fill(180. / TMath::Pi() * evts[i].PhiGen());
            double tmp = evts[i].PhiGen();
            if(tmp > TMath::Pi()) tmp -= TMath::Pi() * 2;
            points -> Fill((180. / TMath::Pi()) * (tmp) * TMath::Sin(HiSC_the), (180. / TMath::Pi()) * (evts[i].ThetaGen()));
        }
    }
    vector<Event> evtsMC;
    ReadTreeOfEvts(FileName2, evtsMC);
    for(int i = 0; i < (int) evtsMC.size(); i++)
    {
        if(evtsMC[i].NStations() >= N_STN_MIN && evtsMC[i].MatchHiSC() == 1){
            thetaDistrMC->Fill(180. / TMath::Pi() * evtsMC[i].ThetaGen());
            phiDistrMC->Fill(180. / TMath::Pi() * evtsMC[i].PhiGen());
        }
    }

    if(draw)
    {
        auto c15 = new TCanvas("c15", "");
        phiDistrMC -> SetLineColor(kRed);
        phiDistrMC -> Scale(1. / phiDistrMC -> GetEntries());
        phiDistrHiSC -> Scale(1. / phiDistrHiSC -> GetEntries());
        phiDistrHiSC -> SetMaximum(phiDistrMC -> GetMaximum() + 0.005);
        phiDistrMC -> Draw("HIST same");
        phiDistrHiSC -> Draw("same");
        auto c16 = new TCanvas("c16", "");
        thetaDistrMC -> SetLineColor(kRed);
        thetaDistrMC -> Scale(1. / thetaDistrMC -> GetEntries());
        thetaDistrHiSC -> Scale(1. / thetaDistrHiSC -> GetEntries());
        thetaDistrMC -> SetMaximum(/*thetaDistrHiSC -> GetMaximum() + 0.005*/ 0.04);
        phiDistrMC -> SetMaximum(0.025);

        thetaDistrMC -> Draw("HIST");
        thetaDistrHiSC -> Draw("same");

        auto c9 = new TCanvas("c9","");
        points -> Draw("CONT4Z");

    }
    return 0;
}