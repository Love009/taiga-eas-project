#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;
#include "TMath.h"

const int N_STN_TOT = 19.;
const char *FILE_TO_READ = "analysis_results/Gr+Hisc010120.root";
const int N_ST_MIN = 3;
double RMS[N_STN_TOT - 1];

int FindSeveralDT(const char *FileName = FILE_TO_READ, int draw = 1)
{
    vector<Event> evts;
    ReadTreeOfEvts(FileName, evts);
    TH1D *DT[N_STN_TOT - 1];
    auto c2 = new TCanvas("c2", "");
    c2 -> Divide(4, 5);
    string name = "DT";
    string title = "Delta t, ns, sts ";
    for(int l = 0; l < N_STN_TOT - 1; l++)
    {
        string s1 = name + to_string(l + 1);
        string s2 = title+ to_string(l+1) + "-" + to_string(l + 2);
        DT[l] = new TH1D(s1.c_str(), s2.c_str(), 200.,  -3000.,  3000.);
       // cout << s1 << "  " << s2 << endl;
        int sti, stj;
        int stnum1 = l + 1;
        int stnum2 = l + 2;
        for(int i = 0; i < (int) evts.size(); i++)
        {
            sti = -1;
            stj = -1;
            for(int j = 0; j < evts[i].NStations(); j++)
            {
                if(evts[i].MuStn(j).ID() == stnum1) sti = j;
                if(evts[i].MuStn(j).ID() == stnum2) stj = j;
            }
            if(sti != -1 && stj != -1 && evts[i].NStations() >= N_ST_MIN /*&& evts[i].MatchHiSC() == 1*/)
            {
                double ts1 = evts[i].MuStn(sti).Trow()[0] * 3600 + evts[i].MuStn(sti).Trow()[1] * 60 + evts[i].MuStn(sti).Trow()[2];
                double t1 = evts[i].MuStn(sti).Trow()[3] + 1000000000. * ts1;
                double ts2 = evts[i].MuStn(stj).Trow()[0] * 3600 + evts[i].MuStn(stj).Trow()[1] * 60 + evts[i].MuStn(stj).Trow()[2];
                double t2 = evts[i].MuStn(stj).Trow()[3] + 1000000000. * ts2;
                DT[l] -> Fill(t2 - t1);
            }
        }
        if(draw)
        {
            gStyle->SetStatW(0.3);
            gStyle->SetStatH(0.4);
            c2 -> cd(l + 1);
            double mean = DT[l] -> GetMean();
            double rms_temp = DT[l] -> GetRMS();
            TF1 *f = new TF1("f", "gaus");
            //if(l != 6) DT[l] -> Fit("f", "Q","", -3000., 3000.);
            /*else*/ DT[l] -> Fit("f", "Q","", mean - 6 * rms_temp, mean + 6 * rms_temp);
            double my_rms = f-> GetParameter(2);
            cout.precision(7);
            cout << "DT_RMS = " << my_rms;
            cout << "\tsts: " << l+1 << "-" << l + 2 << endl;
            RMS[l] = my_rms;
            DT[l] -> GetYaxis()->SetLabelSize(0.06);
            DT[l] -> GetXaxis()->SetLabelSize(0.06);
            gStyle -> SetTitleFontSize(0.09);
            DT[l] -> Draw();
        }
    }
    //delete c2;
    //for(int l = 0; l < N_STN_TOT - 1; l++) delete DT[l];
    return 0;
}

/*int Compare()
{
    FindSeveralDT("analysis_results/Gr+Hisc301219.root", 1);
    FindSeveralDT("analysis_results/Gr+Hisc210220.root", 2);
    cout << "\n\n" << endl;
    for(int u = 0; u < 18; u++) cout << RMS2[u] - RMS1[u] << "\tsts: " << u + 1 << "-" << u + 2 << endl;
    return 0;
}*/


int Write_to_file(const char * FileName, const char * FileToWrite)
{
    ofstream out(FileToWrite);
    out << "stID1,stID2\tRMS2-RMS1, ns" << "\n";
    FindSeveralDT(FileName, 1);
    for(int l = 0; l < N_STN_TOT - 1; l++)
    {
        out << l + 1 << "-" << l + 2 << "\t" << std::setprecision(6) << RMS[l] << endl;
    }
    out.close();
    return 0;
}

int main1()
{
    Write_to_file("analysis_results/Gr+Hisc010120.root", "analysis_results/rms_dt_gr-gr/010120.txt");
    Write_to_file("analysis_results/Gr+Hisc130220.root", "analysis_results/rms_dt_gr-gr/130220.txt");
    Write_to_file("analysis_results/Gr+Hisc150220.root", "analysis_results/rms_dt_gr-gr/150220.txt");
    Write_to_file("analysis_results/Gr+Hisc160220.root", "analysis_results/rms_dt_gr-gr/160220.txt");
    Write_to_file("analysis_results/Gr+Hisc170220.root", "analysis_results/rms_dt_gr-gr/170220.txt");
    Write_to_file("analysis_results/Gr+Hisc200220.root", "analysis_results/rms_dt_gr-gr/200220.txt");
    Write_to_file("analysis_results/Gr+Hisc210220.root", "analysis_results/rms_dt_gr-gr/210220.txt");
    Write_to_file("analysis_results/Gr+Hisc250220.root", "analysis_results/rms_dt_gr-gr/250220.txt");
    Write_to_file("analysis_results/Gr+Hisc270220.root", "analysis_results/rms_dt_gr-gr/270220.txt");

    return 0;
}

/*int Draw_vs_time()//or write to the txt file for each station and date, and after that draw...
{
    TH1D *DTvsTime[18];
    auto c3 = new TCanvas("c3", "");
    c3 -> Divide(4, 5);
    string name = "DTvsTime";
    string title = "Delta t, ns, sts ";
    for(int l = 0; l < 18; l++) {
        string s1 = name + to_string(l + 1);
        string s2 = title + to_string(l + 1) + "-" + to_string(l + 2);
        DT[l] = new TH1D(s1.c_str(), s2.c_str(), 200., -3000., 3000.);
    }
    return 0;
}*/
