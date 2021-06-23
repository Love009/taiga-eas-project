#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;

static const int    N_STN_MIN =   3;//for sistematic errors analysis N_STN_MIN = 2 is enough
const char*         HiSCORE_DATA_ADDRESS   = "HiSCORE_data/2019-2020/out_2611.dat";
const char*         GRANDE_DATA_ADDRESS    = "analysis_results/Gr261119.03.root";
const char*         OUTPUT_FILE_NAME       = "analysis_results/test.root";

long int    DT_WINDOW = 100000;// ns
TH1D *h0   = new TH1D("h0", "dT_Gr_HiSC", 2000,  -1000,  1000);

int DrawADCsignals(int evtNum = 0, const char *grDat = GRANDE_DATA_ADDRESS)
{
    vector<Event> evts, evtsNew;
    ReadTreeOfEvts(grDat, evts);
    TCanvas *c = new TCanvas("c", "Signal from channel");
     c->Divide(4, 3);
     for(int u = 0; u < 8; u++)//8 channels
     {
         TH1I *sig = new TH1I("sig", "Signal from channel", evts[evtNum].MuStn(0).ADC(u).size(), 0,
                              evts[evtNum].MuStn(0).ADC(u).size());
         sig->GetYaxis()->SetRangeUser(2000, 2200);
         for (int yy = 0; yy < (int) evts[evtNum].MuStn(0).ADC(u).size(); yy++) sig->SetBinContent((yy + 1),
                                                                                              evts[evtNum].MuStn(0).ADC(
                                                                                                      u)[yy]);
         c->cd(u + 1);
         sig->Draw("SAME");
     }
    return 0;
}
