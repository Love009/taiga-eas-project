#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;

static const int    N_STN_TOT =  19;
static const double HiSC_CALIBR =-1537.;/*136758000.+230.;*/
const char*         FILE_NAME   = "analysis_results/Gr+Hisc241119.root";
TH1D *h0   = new TH1D("h0", "dT_Gr_HiSC, ns", 1000,  (-10000 /*- HiSC_CALIBR*/),  (10000 /*- HiSC_CALIBR*/));
TH1D *h1   = new TH1D("h1", "dT_Gr_HiSC, ns", 1000,  (-10000 /*- HiSC_CALIBR*/),  (10000 /*- HiSC_CALIBR*/));
TH2D *corrx   = new TH2D("corrx", "correlation (dtGr-Hisc, x)", 1200, -600, 600, 200,  (-10000 - HiSC_CALIBR)/1000000.,  (10000 - HiSC_CALIBR)/1000000.);
TH2D *corry   = new TH2D("corry", "correlation (dtGr-Hisc, y)", 1200, -600, 600, 200,  (-10000 - HiSC_CALIBR)/1000000.,  (10000 - HiSC_CALIBR)/1000000.);
TH2D *corrr   = new TH2D("corrr", "correlation (dtGr-Hisc, r)", 1600, -800, 800, 200,  (-10000 - HiSC_CALIBR)/1000000.,  (10000 - HiSC_CALIBR)/1000000.);
TH2D *corrthe   = new TH2D("corrthe", "correlation (dtGr-Hisc, theta)", 90, 0, 90, 200,  (-10000 - HiSC_CALIBR)/1000000.,  (10000 - HiSC_CALIBR)/1000000.);
TH2D *corrphi   = new TH2D("corrphi", "correlation (dtGr-Hisc, phi)", 180, 0, 180, 200,  (-10000 - HiSC_CALIBR)/1000000.,  (10000 - HiSC_CALIBR)/1000000.);
//===========================================================================================

long int GetDT(Event &evt)
{
    long int dT0 = (evt.EvtTime()[0] - evt.EvtTimeHiSC()[0]) * 3600 +(evt.EvtTime()[1] - evt.EvtTimeHiSC()[1]) * 60 + (evt.EvtTime()[2] - evt.EvtTimeHiSC()[2]);
    long int dT = dT0 * 1000000000 + evt.EvtTime()[3] - evt.EvtTimeHiSC()[3];
    return dT;
}

int DrawDTGrHisc(const char* fileName = FILE_NAME, int hiscCalibr = HiSC_CALIBR)
{
    vector<Event> evts;
    ReadTreeOfEvts(fileName, evts);
    for(int k = 0; k < (int) evts.size();k++)
    {
        int here = 0, ten=0, nine = 0;
        for(int u = 0; u < evts[k].NStations(); u++) if(evts[k].MuStn(u).ID() == 10) {here = 1;ten=u;}
        if(here==1)for(int u = 0; u < evts[k].NStations(); u++)if(u!=ten)
        {
            double dr = sqrt(pow((evts[k].MuStn(u).X() - evts[k].MuStn(ten).X()),2) + pow((evts[k].MuStn(u).Y() - evts[k].MuStn(ten).Y()),2));
           // if(dr > 250) here = 0;
            if(evts[k].MuStn(u).ID() != 9) {here =0;nine=u;}
        }
        if(here == 1)
        {
            h1->Fill(evts[k].MuStn(ten).T()-evts[k].MuStn(nine).T());
        long int dT = GetDT(evts[k]);
	cout << dT + HiSC_CALIBR<< endl;
        h0->Fill(dT+HiSC_CALIBR);
        corrx->Fill(evts[k].MuStn(0).X(), (dT - HiSC_CALIBR)/1000000.);
        corry->Fill(evts[k].MuStn(0).Y(), (dT - HiSC_CALIBR)/1000000.);
        corrr->Fill(sqrt(pow(evts[k].MuStn(0).X(),2)+pow(evts[k].MuStn(0).Y(),2)), (dT - HiSC_CALIBR)/1000000.);
        corrthe->Fill(evts[k].ThetaGen() * 57, (dT - HiSC_CALIBR)/1000000.);
        corrphi->Fill(evts[k].PhiGen() * 57, (dT - HiSC_CALIBR)/1000000.);
        }
    }
    auto c1 = new TCanvas("c1", "");
    h0->Draw();
    
    auto c8 = new TCanvas("c8", "");
    h1->Draw();
    
   // auto c2 = new TCanvas("c2", "");
 //   corrr->Draw("BOX");
    //auto c3 = new TCanvas("c3", "");
   // corrx->Draw("BOX");
    //auto c4 = new TCanvas("c4", "");
   // corry->Draw("BOX");
    auto c5 = new TCanvas("c5", "");
    corrthe->Draw();
    auto c6 = new TCanvas("c6", "");
    corrphi->Draw();
    return 0;
}
