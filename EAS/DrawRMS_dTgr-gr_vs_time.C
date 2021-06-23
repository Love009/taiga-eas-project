#include "common.h"
#include "ReadTreeOfEvts.h"
#include <iostream>
#include <dirent.h>
#include <sys/types.h>
using namespace std;

const int TIME_WINDOW = 5; // *5 ns
//vector<double> MeanCH0;

static const int    N_STN_TOT =  19;
double              MY_RMS[N_STN_TOT - 1];
//===============================================================================================//
//Initialize station channel pedestals
void Initialize(const char* fileName)
{
    double rms_;
    string  title, rest;
    int i = 0;
    ifstream in(fileName);
    getline(in, title);
    while (1) {
        in >> rest >> rms_;
        if (!in.good()) break;
        MY_RMS[i] = rms_;
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
        dates[i][0] = stoi(fileNames[i].substr(0,2));
        dates[i][1] = stoi(fileNames[i].substr(2,2));
        dates[i][2] = stoi(fileNames[i].substr(4,2));
    }
    PreSort(dates, n, fileNames, 2);
    PreSort(dates, n, fileNames, 1);
    PreSort(dates, n, fileNames, 0);
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CalcMean(TH1D *hist, int n)
{
    double mean = 0;
    int counter = 0;
    for(int i = 1; i <= n; i++){
        double tmp0 = hist -> GetBinContent(i);
        if(tmp0 != 0){ mean += tmp0; counter++;}
    }
    mean = mean / counter;
    return mean;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int DrawRMS_vs_time(int sepDrStNum = 1)
{
    TCanvas *c2 = new TCanvas("c2", "RMS");
    c2 -> Divide(4, 5);
    TH1D *my_rms[N_STN_TOT - 1];
    string name = "my_rms";
    string title = "RMS(delta t), ns, sts ";

    vector<string> fileNames;
    list_dir("analysis_results/rms_dt_gr-gr", fileNames);
    int n = (int)fileNames.size();
    for(int i = 0; i < n; i++) if(fileNames[i][0] == '.') {fileNames.erase(fileNames.begin() + i); i--; n--;}
    SortFiles(fileNames);
    n = (int)fileNames.size();
    for(int l = 0; l < N_STN_TOT - 1; l++) {
        string s1 = name + to_string(l + 1);
        string s2 = title + to_string(l + 1) + "-" + to_string(l + 2);
        my_rms[l] = new TH1D(s1.c_str(), s2.c_str(), n, 0, n);
    }

    for(int i = 0; i < n; i++)
    {
        cout << fileNames[i] << endl;
        Initialize(("analysis_results/rms_dt_gr-gr/"+fileNames[i]).c_str());
        for(int stID = 1; stID < N_STN_TOT; stID++) my_rms[stID - 1] -> SetBinContent(i + 1, MY_RMS[stID - 1]);
    }

    for(int stID = 1; stID < N_STN_TOT; stID++) {
        c2->cd(stID);
        my_rms[stID - 1] -> SetMarkerStyle(7);
        my_rms[stID - 1] -> SetMarkerSize(2);
        my_rms[stID - 1] -> SetMarkerColor(kBlue);
        my_rms[stID - 1] -> GetYaxis()->SetLabelSize(0.06);
        my_rms[stID - 1] -> GetXaxis()->SetLabelSize(0.06);
        int n1 = my_rms[stID - 1]->GetXaxis()->GetNbins();
        double mean = CalcMean(my_rms[stID - 1],n1);
        my_rms[stID - 1] -> SetAxisRange(mean - 300., mean + 300.,"Y");
        gStyle -> SetTitleFontSize(0.09);
        my_rms[stID - 1] -> SetStats(0);
        my_rms[stID - 1] -> Draw("PL");
    }

    return 0;
}

