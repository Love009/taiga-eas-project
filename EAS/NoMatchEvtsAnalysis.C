//prints efficiency of HiSCORE
#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;

const char *FILE_TO_ANALYZE = "analysis_results/for_EAS_center/Gr+Hisc170220Recon.root";
const char*         HiSCORE_DATA_ADDRESS   = "analysis_results/out_1702.dat";
int T1[3], T2[3];

static const int    N_STN_MIN =   3;
static const int    N_STN_TOT =  19;
TH1D *nTrigSts = new TH1D("nTrigSts",   "Number of stations triggered in event, matchHISC=1", (N_STN_TOT ),  1,  N_STN_TOT + 1);
TH1D *nTrigSts2 = new TH1D("nTrigSts2",   "Number of stations triggered in event, matchHISC=0", (N_STN_TOT ),  1,  N_STN_TOT + 1);
TH1D *phihi = new TH1D("phihi",   "phi NO match", 360./2.,  -360,  360);
TH1D *theta = new TH1D("theta",   "theta NO match", 180./2.,  -180,  180);

TH1D *phihiMatched = new TH1D("phihiMatched",   "phi Matched", 360./2.,  -360,  360);
TH1D *thetaMatched = new TH1D("thetaMatched",   "theta Matched", 180./2.,  -180,  180);
//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
int HiscTime(const char* hiscDat)
{
    vector<Event> evts, evtsNew;
    string line, lineLast;
    ifstream hiscore(hiscDat);
    if (! hiscore.is_open()) return -1;
    int count = 0;
    const char u = ':';
    while(getline(hiscore, line))
    {
        count++;
        int k = 0;
        while(char(line[k]) != u)  k++;
        if(count == 1)
        {
            T1[0] = stoi(line.substr(k - 2, 9));
            T1[1] = stoi(line.substr(k + 1, 2));
            T1[2] = stoi(line.substr(k + 4, 2)) + 1;
        }
        lineLast = line;
    }
    int k = 0;
    while(char(lineLast[k]) != u)  k++;
    T2[0] = stoi(lineLast.substr(k - 2, 9));
    T2[1] = stoi(lineLast.substr(k + 1, 2));
    T2[2] = stoi(lineLast.substr(k + 4, 2));
    hiscore.close();
    cout << "HiSCORE time interval: " << T1[0] << ":" << T1[1] <<  ":" << T1[2] <<" - " << T2[0] << ":" << T2[1] <<  ":" << T2[2] << endl;
    return 0;
}
//------------------------------------------------------------------------------------------------------
int CheckCenter600(Event &evt)
{
    if(sqrt(evt.EAScentersGen()[0] * evt.EAScentersRec()[0] + evt.EAScentersRec()[1] * evt.EAScentersRec()[1]) <= 400.) return 1;
    else return 0;
}
//------------------------------------------------------------------------------------------------------
int main1(const char *fileToAnalyze = FILE_TO_ANALYZE, const char *hiscDat = HiSCORE_DATA_ADDRESS)
{
    nTrigSts->Reset();
    nTrigSts2->Reset();
    HiscTime(hiscDat);
    vector<Event> evts;
    int n3sts_tot = 0, n3sts = 0;
    int t1 = T1[0] * 60 * 60 + T1[1] * 60 + T1[2];
    int t2 = T2[0] * 60 * 60 + T2[1] * 60 + T2[2];
    double the_hisc_decl = 25. * TMath::Pi() / 180.;
    double hisc_res_size = 25. * TMath::Pi() / 180.;
    ReadTreeOfEvts(fileToAnalyze, evts);
    for(int i = 0; i < (int)evts.size(); i++)
    {
        double the_rec = evts[i].ThetaRec();
        double phi_rec = evts[i].PhiRec();
        double temp = TMath::ACos(sin(the_hisc_decl) * sin(the_rec) * cos(phi_rec) + cos(the_hisc_decl) * cos(the_rec));
        int flag600 = CheckCenter600(evts[i]);
        if(evts[i].MatchHiSC() == 1 && evts[i].NStations() >= N_STN_MIN) {
            long int t = evts[i].MuStn(0).Trow()[0] * 60 * 60 + evts[i].MuStn(0).Trow()[1] * 60 + evts[i].MuStn(0).Trow()[2];
            double the = evts[i].ThetaRec();
            if (t >= t1 && t <= t2 && evts[i].Chi2() < 10 && evts[i].PhiRec() != -100. && evts[i].ThetaRec() != -100. && flag600 == 1)
            {
                if(temp <= hisc_res_size)
                {
                    nTrigSts->Fill(evts[i].NStations());
                    phihiMatched->Fill(57.3 * evts[i].PhiRec());
                    thetaMatched->Fill(57.3 * evts[i].ThetaRec());
                }
                else cout << "\t\t\t" << temp << endl;
            }
           // else cout << t << "(" << t2 << "-" << t1 << ")\n";

        }
        else if(evts[i].MatchHiSC() == 0 && evts[i].NStations() >= N_STN_MIN)
        {
            long int t = evts[i].MuStn(0).Trow()[0] * 60 * 60 + evts[i].MuStn(0).Trow()[1] * 60 + evts[i].MuStn(0).Trow()[2];
            if(t >= t1 && t <= t2 && evts[i].Chi2() < 10 && evts[i].PhiRec() != -100. && evts[i].ThetaRec() != -100. && flag600 == 1)
            {
                if(temp <= hisc_res_size)
                {
                    nTrigSts2->Fill(evts[i].NStations());
                    phihi->Fill(57.3 * evts[i].PhiRec());
                    theta->Fill(57.3 * evts[i].ThetaRec());
                }
            }
        }
        long int t = evts[i].MuStn(0).Trow()[0] * 60 * 60 + evts[i].MuStn(0).Trow()[1] * 60 + evts[i].MuStn(0).Trow()[2];
        if(evts[i].NStations() == 3 && t >= t1 && t <= t2) n3sts_tot++;
        if(evts[i].NStations() >= 3 && evts[i].PhiRec() != -100. && evts[i].ThetaRec() != -100. && t >= t1 && t <= t2 && evts[i].Chi2() < 10 && temp <= hisc_res_size * hisc_res_size) n3sts++;


    }


    int n1, n2, eff[2];
    n1 = nTrigSts -> GetBinContent(3);
    n2 = nTrigSts2 -> GetBinContent(3);
    eff[0] = int(100. * n1 / (n1 + n2));
    n1 = 0;
    n2 = 0;
    for(int u = 4; u <= 19; u++)
    {
        n1 += nTrigSts -> GetBinContent(u);
        n2 += nTrigSts2 -> GetBinContent(u);
    }

    eff[1] = int(100. * n1 / (n1 + n2));
    cout << eff[0] << "\n" << eff[1] << endl;
   /* for(int j = 3; j <= 6; j++)
    {
        int n1, n2, eff;
        n1 = nTrigSts -> GetBinContent(j);
        n2 = nTrigSts2 -> GetBinContent(j);
        eff = 100. * n1 / (n1 + n2);
        cout << eff << endl;
    }*/
    auto c3 = new TCanvas("c3", "");
    c3 -> Divide(2, 1);
    c3 -> cd(1);
    phihiMatched -> Draw();
    c3 -> cd(2);
    phihi -> Draw();
    auto c4 = new TCanvas("c4", "");
    c4 -> Divide(2, 1);
    c4 -> cd(1);
    thetaMatched -> Draw();
    c4 -> cd(2);
    theta -> Draw();


    auto c2 = new TCanvas("c2", "");
    c2 -> Divide(2, 1);
    c2 -> cd(1);
    nTrigSts -> Draw();
    c2 -> cd(2);
    nTrigSts2 -> SetMaximum(nTrigSts -> GetMaximum() + 50.);
    nTrigSts2 -> Draw();
    cout << "3 sts evts num tot = " << n3sts_tot << endl;
    cout << "3 sts evts num selected = " << n3sts << endl;
    return 0;
}


