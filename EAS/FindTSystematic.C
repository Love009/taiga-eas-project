//systematic shifts between Grande stations - prototype for func in the file Recon_new.C
#include "common.h"
#include "ReadTreeOfEvts.h"
#include <iostream>
#include <dirent.h>
#include <sys/types.h>
using namespace std;

static const double c = 0.3; //speed of light, m/ns
static const int    N_STN_TOT =  19;
static const int    N_STN_MIN = 2;
vector<double> MuonX;
vector<double> MuonY;
vector<double> MuonZ;
const char*    FILE_NAME    = "ToyMC.root";
double statist;
static const int    MATCH_HiSCORE = 1;
static const int    N_STN_MIN_SYS = 2;
//===============================================================================================//
//Initualize stations coordinates and time errors (vectors MuonX, Y, Z, Terr)
// par = 1 for onground stations; par = 2 for underground
void Initialize(int par = 1)
{
    double x, y, z, sigma;
    string  title, rest;
    int i = 1;
    ifstream in("coordinatesWithTimeErr.txt");
    getline(in, title);
    while (1) {
        in >> rest >> x >> y >> z >> sigma;// sigma = t_error
        if (!in.good()) break;
        if(i % 2 == par % 2){//2%2 for getting underground counter coordinstes, 1%2 for onground
            MuonX.push_back(x);
            MuonY.push_back(y);
            MuonZ.push_back(z);
        }
        i++;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double tSyst[N_STN_TOT];
double *FindTSystematic(const char *fileName = FILE_NAME)
{
    TMatrixD e(N_STN_TOT,1);
    for(int u = 0; u < N_STN_TOT; u++) e(u,0) = 1.;
    int n_useful = 0, n_selected = 0;
    Initialize(1);// par = 1 <- onground stations coordinates
    TMatrixD R(N_STN_TOT, 3);
    for(int i = 0; i < N_STN_TOT; i++)
    {
        R(i,0) = MuonX[i];
        R(i,1) = MuonY[i];
        R(i,2) = MuonZ[i];
    }

    vector<Event> evts;
    ReadTreeOfEvts(fileName, evts);
    TMatrixD dt(N_STN_TOT,1);
    TMatrixD res(N_STN_TOT,1);
    TMatrixD M_tilda_tot(N_STN_TOT, N_STN_TOT);
    for(int u = 1; u < N_STN_TOT; u++)for(int uu = 1; uu < N_STN_TOT; uu++) M_tilda_tot(u,uu) = 0.;
    for(int u = 1; u < N_STN_TOT; u++) res(u,0) = 0.;
    for(int k = 0; k < (int)evts.size(); k++) {
        if (evts[k].MatchHiSC() == MATCH_HiSCORE)
        {
            TMatrixD M_k(N_STN_TOT, N_STN_TOT);
            for (int jj = 0; jj < N_STN_TOT; jj++) for (int ii = 0; ii < N_STN_TOT; ii++) M_k(ii, jj) = 0.;
            TMatrixD M_tilda_k(N_STN_TOT, N_STN_TOT);
            TMatrixD n_k(1, 3);
            TMatrixD t_k(N_STN_TOT, 1);
            TMatrixD n_k_R(N_STN_TOT, 1);
            int mk = 0;

            if (evts[k].NStations() == N_STN_MIN_SYS) {
                n_useful++;
                Event evt = evts[k];
                n_k(0, 0) = cos(evt.PhiGen()) * sin(evt.ThetaGen());
                n_k(0, 1) = sin(evt.PhiGen()) * sin(evt.ThetaGen());
                n_k(0, 2) = cos(evt.ThetaGen());
                for (int i = 0; i < N_STN_TOT; i++)
                    n_k_R(i, 0) = n_k(0, 0) * R(i, 0) + n_k(0, 1) * R(i, 1) + n_k(0, 2) * R(i, 2);
                mk = evt.NStations();

                double dtime;
                if (evts[k].NStations() == N_STN_MIN_SYS) {
                    dtime = (evt.MuStn(0).T() - evt.MuStn(1).T());
                    double nr = 0.;
                    for (int ll = 0; ll < 3; ll++)
                        nr += n_k(0, ll) * (R(evt.MuStn(0).ID() - 1, ll) - R(evt.MuStn(1).ID() - 1, ll));
                    dtime -= nr / 0.3;
                }

                for (int j = 0; j < mk; j++) {
                    int stID = evt.MuStn(j).ID();
                    t_k(stID - 1, 0) = evt.MuStn(j).T();//?
                    M_k(stID - 1, stID - 1) = 1.;

                }

                TMatrixD delta(N_STN_TOT, 1);
                delta = M_k * e;
                for (int jj = 0; jj < N_STN_TOT; jj++)
                    for (int ii = 0; ii < N_STN_TOT; ii++) {
                        M_tilda_k(ii, jj) = delta(ii, 0) * delta(jj, 0) / mk;
                    }
                if (abs(dtime) < 500) {
                    n_selected++;
                    M_tilda_tot += M_k - M_tilda_k;
                    res += (M_k - M_tilda_k) * (-1 / c * n_k_R + t_k);
                }
            }
        }
    }

    for(int u = 0; u < N_STN_TOT; u++) M_tilda_tot(0, u) = 1.;
    res(0, 0) = 0.;
    
    //ref point at the station ID=0:
    //for(int u = 0; u < N_STN_TOT; u++) M_tilda_tot(0, u) = 0.;
    // M_tilda_tot(0, 18) = 1.;
    
    
    TMatrixD M_tilda_tot_inv = M_tilda_tot.Invert();
    dt = M_tilda_tot_inv * res;


    double sum = 0., sumabs = 0.;
    for(int i = 0; i < N_STN_TOT; i++) {sum += dt(i,0); sumabs += abs(dt(i,0)); tSyst[i] = dt(i,0);}
    dt.Print();//ns
    // cout << "sum = "<< sum << "\nsumabs = " << sumabs << "\n" << n_useful << " useful events" << endl;
    cout << n_selected << " selected events" << endl;
    cout << n_useful << " useful evts" << endl;
    return tSyst;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Reads file names from the directory path and fills fileNames array
void list_dir(const char *path, vector<string> &fileNames) {
    struct dirent *entry;
    DIR *dir = opendir(path);

    if (dir == NULL) {
        return;
    }
    while ((entry = readdir(dir)) != NULL) fileNames.push_back(entry->d_name);
    closedir(dir);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Swap(string *xp, string *yp)
{
    string temp = *xp;
    *xp = *yp;
    *yp = temp;
}
//------------------------------------------------------------------------//pedestal121019.txt  Gr+Hisc121019.root
void PreSort(int (*dates)[3], int n, vector<string> &fileNames, int part)
{
    for(int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n-i-1; j++)
        {
            if(part == 2)
            {
                if (dates[j][part] > dates[j + 1][part]) {
                    Swap(&fileNames[j], &fileNames[j + 1]);
                    int temp[3];
                    for(int k = 0; k < 3; k++) temp[k] = dates[j][k];
                    for(int k = 0; k < 3; k++) dates[j][k] = dates[j + 1][k];
                    for(int k = 0; k < 3; k++) dates[j + 1][k] = temp[k];
                }
            }
            if(part == 1)
            {
                if (dates[j][part] > dates[j + 1][part] && dates[j][2] == dates[j + 1][2]) {
                    Swap(&fileNames[j], &fileNames[j + 1]);
                    int temp[3];
                    for(int k = 0; k < 3; k++) temp[k] = dates[j][k];
                    for(int k = 0; k < 3; k++) dates[j][k] = dates[j + 1][k];
                    for(int k = 0; k < 3; k++) dates[j + 1][k] = temp[k];
                }
            }
            if(part == 0)
            {
                if (dates[j][part] > dates[j + 1][part] && dates[j][2] == dates[j + 1][2] &&  dates[j][1] == dates[j + 1][1]) {
                    Swap(&fileNames[j], &fileNames[j + 1]);
                    //Swap(dates[j], dates[j + 1]);
                    int temp[3];
                    for(int k = 0; k < 3; k++) temp[k] = dates[j][k];
                    for(int k = 0; k < 3; k++) dates[j][k] = dates[j + 1][k];
                    for(int k = 0; k < 3; k++) dates[j + 1][k] = temp[k];
                }
            }
        }
    }
}
//------------------------------------------------------------------------//
int SortFiles(vector<string> &fileNames)
{
    int n = (int)fileNames.size();
    int dates[n][3];
    for(int i = 0; i < n; i++)
    {
        dates[i][0] = stoi(fileNames[i].substr(7,2));
        dates[i][1] = stoi(fileNames[i].substr(9,2));
        dates[i][2] = stoi(fileNames[i].substr(11,2));
    }
    PreSort(dates, n, fileNames, 2);
    PreSort(dates, n, fileNames, 1);
    PreSort(dates, n, fileNames, 0);
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int DrawTsystematic_vs_time(int sepDrStNum = 1)
{
    vector<string> fileNames;
    list_dir("analysis_results/no_syn2/", fileNames);
    int n = (int)fileNames.size();
    cout << "begin " << n << endl;
    for(int i = 0; i < n; i++) if(fileNames[i][0] == '.') {fileNames.erase(fileNames.begin() + i); i--; n--;}
    SortFiles(fileNames);
    n = (int)fileNames.size();
    for(int i = 0; i < n; i++) cout << fileNames[i] << endl;
    
    vector<double> tSyst[N_STN_TOT];
    double statist_[n];
    statist = 0.;
    for(int i = 0; i < n; i++)
    {
        double *temp;
        temp = FindTSystematic(("analysis_results/no_syn2/"+fileNames[i]).c_str());
        statist_[i] = statist;
        for(int j = 0; j < N_STN_TOT; j++) tSyst[j].push_back(*(temp +j));
    }

    TCanvas *c1 = new TCanvas("c1", "");
    //c1 -> Divide(4, 5);
    c1 -> Divide(3, 4);
    for(int stID = 0; stID < 12/*N_STN_TOT*/; stID++)
    {
        c1 -> cd(stID + 1);
        TH1D *dT = new TH1D("dT", "dt systematic, ns", n, 0, n);
        for(int i = 0; i < (int)tSyst[stID].size(); i++){dT -> SetBinContent(i + 1, tSyst[stID][i]/* - tSyst[0][i]*/); /*dT->SetBinError(i + 1, 1 * pow(20000. / statist_[i], 0.5));*/}
        dT -> SetMarkerStyle(7);
        dT -> SetMarkerSize(2);
        dT -> SetMarkerColor(kBlue);
        dT -> SetAxisRange(-300., 300.,"Y");
       // dT -> GetXaxis() -> SetTitle("day");
      //  dT -> GetYaxis() -> SetTitle("dt syst");
        gStyle->SetOptStat(0);

        dT -> GetYaxis()->SetLabelSize(0.08);
        //dT -> GetYaxis()->SetTitleSize(0.03);
        dT -> GetXaxis()->SetLabelSize(0.08);
       // dT -> GetXaxis()->SetTitleSize(0.03);
        gStyle -> SetTitleFontSize(0.10);
        dT -> Draw("PL SAME E");
    }
    return 0;
    
}
