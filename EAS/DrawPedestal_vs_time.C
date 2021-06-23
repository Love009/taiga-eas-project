#include "common.h"
#include "ReadTreeOfEvts.h"
#include <iostream>
#include <dirent.h>
#include <sys/types.h>
using namespace std;

const int TIME_WINDOW = 5; // *5 ns
vector<double> MeanCH0;
vector<double> MeanCH2;
vector<double> SigmaCH0;
vector<double> SigmaCH2;
TH1D *signalStart = new TH1D("signalStart","signalStart", 200,  0,  200); //hist for linear fit
const char *INITIALIZATION_FILE_NAME = "analysis_results/pedestals/pedestal130220.txt";
static const int    N_STN_TOT =  19;
//===============================================================================================//
//Initialize station channel pedestals
void Initialize(const char* fileName = INITIALIZATION_FILE_NAME)
{
    double meanCH0, meanCH2, sigmaCH0, sigmaCH2;
    MeanCH0.clear();
    MeanCH2.clear();
    SigmaCH0.clear();
    SigmaCH2.clear();
    string  title, rest;
    int i = 1;
    ifstream in(fileName);
    getline(in, title);
    while (1) {
        in >> rest >> meanCH0 >> sigmaCH0 >> meanCH2 >> sigmaCH2;
        if (!in.good()) break;
        MeanCH0.push_back(meanCH0);
        MeanCH2.push_back(meanCH2);
        SigmaCH0.push_back(sigmaCH0);
        SigmaCH2.push_back(sigmaCH2);
        i++;
    }
}
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
//------------------------------------------------------------------------//
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
                    //Swap(dates[j], dates[j + 1]);
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
                    //Swap(dates[j], dates[j + 1]);
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
        dates[i][0] = stoi(fileNames[i].substr(8,2));
        dates[i][1] = stoi(fileNames[i].substr(10,2));
        dates[i][2] = stoi(fileNames[i].substr(12,2));
    }
    PreSort(dates, n, fileNames, 2);
    PreSort(dates, n, fileNames, 1);
    PreSort(dates, n, fileNames, 0);
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int DrawPedestal_vs_time(int sepDrStNum = 1)
{   
    TCanvas *c1 = new TCanvas("c1", "Pedestal");
    TCanvas *c2 = new TCanvas("c2", "RMS");
    c1 -> Divide(4, 5);
    c2 -> Divide(4, 5);
    TCanvas *can1 = new TCanvas("can1", "Pedestal  ");
    TCanvas *can2 = new TCanvas("can2", "RMS  ");
    vector<string> fileNames;
    list_dir("analysis_results/pedestals", fileNames);
    int n = (int)fileNames.size();
    for(int i = 0; i < n; i++) if(fileNames[i][0] == '.') {fileNames.erase(fileNames.begin() + i); i--; n--;}
    
    SortFiles(fileNames);
    n = (int)fileNames.size();
    
    for(int stID = 1; stID < 20; stID++)
    {
        TH1D * pedestal0 = new TH1D("pedestal0", "Pedestal", n, 0, n);
        TH1D * sigma0 = new TH1D("sigma0", "Sigma", n, 0, n);
        TH1D * pedestal2 = new TH1D("pedestal2", "Pedestal", n, 0, n);
        TH1D * sigma2 = new TH1D("sigma2", "Sigma", n, 0, n);
        for(int i = 0; i < n; i++)
        {
            //cout << fileNames[i] << endl;
            Initialize(("analysis_results/pedestals/"+fileNames[i]).c_str());
            double mean = MeanCH0[stID - 1], rms = SigmaCH0[stID - 1];
            pedestal0 -> SetBinContent(i + 1, mean);
            sigma0 -> SetBinContent(i + 1, rms);
            mean = MeanCH2[stID - 1];
            rms = SigmaCH2[stID - 1];
            pedestal2 -> SetBinContent(i + 1, mean);
            sigma2 -> SetBinContent(i + 1, rms);
        }
        c1 -> cd(stID);
        pedestal0 -> SetMarkerStyle(7);
        pedestal0 -> SetMarkerSize(2);
        pedestal0 -> SetMarkerColor(kBlue);
        pedestal0 -> SetAxisRange(2038., 2058.,"Y");
        gStyle->SetOptStat(0);
        pedestal0 -> Draw("PL");
        pedestal2 -> SetMarkerStyle(7);
        pedestal2 -> SetMarkerSize(2);
        pedestal2 -> SetMarkerColor(kRed);
        pedestal2 -> SetLineColor(kRed);
        pedestal2 -> SetAxisRange(2038., 2058.,"Y");
        pedestal2 -> Draw("PL SAME");
        if(stID == sepDrStNum)
        {
            can1 -> cd();
            //pedestal0 -> SetLineColor(56 + 4 * stID);
            pedestal0 -> GetXaxis()->SetTitle("day");
            pedestal0 -> GetYaxis()->SetTitle("mean, counts");
            pedestal0 -> SetMarkerStyle(20);
            pedestal0 -> SetMarkerSize(1);
            pedestal0 -> SetMarkerColor(kBlue);
            pedestal0 -> Draw("PL SAME");
            pedestal2 -> GetXaxis()->SetTitle("day");
            pedestal2 -> GetYaxis()->SetTitle("mean, counts");
            pedestal2 -> SetMarkerStyle(20);
            pedestal2 -> SetMarkerSize(1);
            pedestal2 -> SetMarkerColor(kRed);
            pedestal2 -> Draw("PL SAME");
        }
        
        c2 -> cd(stID);
        sigma0 -> SetMarkerStyle(7);
        sigma0 -> SetMarkerSize(2);
        sigma0 -> SetMarkerColor(kBlue);
        sigma0 -> SetAxisRange(-4., 4.,"Y");
        sigma0 -> Draw("PL");
        sigma2 -> SetMarkerStyle(7);
        sigma2 -> SetMarkerSize(2);
        sigma2 -> SetMarkerColor(kRed);
        sigma2 -> SetLineColor(kRed);
        sigma2 -> SetAxisRange(-10., 10.,"Y");
        sigma2 -> Draw("PL SAME");
        if(stID == sepDrStNum)
        {
            can2 -> cd();
            sigma0 -> GetXaxis()->SetTitle("day");
            sigma0 -> GetYaxis()->SetTitle("RMS, counts");
            sigma0 -> SetMarkerStyle(20);
            sigma0 -> SetMarkerSize(1);
            sigma0 -> SetMarkerColor(kBlue);
            sigma0 -> Draw("PL SAME");
            sigma2 -> GetXaxis()->SetTitle("day");
            sigma2 -> GetYaxis()->SetTitle("RMS, counts");
            sigma2 -> SetMarkerStyle(20);
            sigma2 -> SetMarkerSize(1);
            sigma2 -> SetMarkerColor(kRed);
            sigma2 -> Draw("PL SAME");
        }
    }
    return 0;
}

