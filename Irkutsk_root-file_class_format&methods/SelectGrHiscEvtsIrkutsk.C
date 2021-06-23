#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;

static const int    N_STN_MIN =   3;//for sistematic errors analysis N_STN_MIN = 2 is enough
//static const int    HiSC_CALIBR = 0;//minutes
const char*         HiSCORE_DATA_ADDRESS   = "HiSCORE_data/2019-2020/out_2611.dat";
const char*         GRANDE_DATA_ADDRESS    = "analysis_results/Gr261119.root";
const char*         OUTPUT_FILE_NAME       = "analysis_results/test.root";

long int    DT_WINDOW = 10000;// ns
TH1D *h0   = new TH1D("h0", "dT_Gr_HiSC", 2000,  -1000,  1000);
//long int    DT_WINDOW = 7000000000; // ns
//TH1D *h0   = new TH1D("h0", "dT_Gr_HiSC", 14000,  -7000,  7000);
//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int WriteToFile(vector<Event> &evtsNew, const char* fileToWrite)
{
    TFile f(fileToWrite,"recreate");
    TTree Muon("Muon","");
    gROOT->ProcessLine("#include <vector>");
    Event *event = new Event();
    Muon.Branch("events", &event);
    Muon.SetAutoSave(1000000);
    for(int i = 0; i < (int)evtsNew.size(); i++)
    {
        event = & evtsNew[i];
        Muon.Fill();
    }
    f.Write();
    f.Close();
    cout<< (int)evtsNew.size() << endl;
    return 0;
}
int WriteEAScenter(string &line, Event &evt, int count)
{
    const char u2 = '.';
    double centerX, centerY;
    int flag = 0;
    int k = count;
    do{
        k++;
        if (char(line[k]) == u2) flag++;
    }while(flag < 6);
    if (line[k - 2] != ' ' && line[k - 3] != ' ' && line[k - 4] != ' ')
        centerX = stod(line.substr(k - 4, 4));
    else if (line[k - 2] != ' ' && line[k - 3] != ' ') centerX = stod(line.substr(k - 3, 3));
    else if (line[k - 2] != ' ') centerX = stod(line.substr(k - 2, 2));
    else centerX = stod(line.substr(k - 1, 1));

    k += 6;
    if (line[k - 2] != ' ' && line[k - 3] != ' ' && line[k - 4] != ' ')
        centerY = stod(line.substr(k - 4, 4));
    else if (line[k - 2] != ' ' && line[k - 3] != ' ') centerY = stod(line.substr(k - 3, 3));
    else if (line[k - 2] != ' ') centerY = stod(line.substr(k - 2, 2));
    else centerY = stod(line.substr(k - 1, 1));
    double coordArr[2] = {centerX, centerY};
    evt.SetXYHiSC(coordArr);
    return 0;
}
int WriteEASdirection(string &line, int k, Event &evt)
{
    evt.SetThetaHiSC(stod(line.substr(k + 20, 5)) / 57.3);
    double pphi = stod(line.substr(k + 26, 6));
    evt.SetPhiHiSC(2 * TMath::Pi() - pphi / 57.3);

    evt.SetRaHiSC(stod(line.substr(k + 89, 6)));
    evt.SetDecHiSC(stod(line.substr(k + 97 , 5)));
    //cout << stod(line.substr(k + 20, 5)) << "  " << evt.RaHiSC() << "  " << evt.DecHiSC() << endl;
    return 0;
}
int GetHiSCtime(const char sep, string line, int &k, int *THiSC)
{
    while(char(line[k]) != sep)  k++;
    long int ns = stol(line.substr(k + 7, 3)) * 1000000 + stol(line.substr(k + 11, 3)) * 1000 + stol(line.substr(k + 15, 3));
    string time = line.substr(k - 2, 9) + to_string(ns);
    THiSC[0] = stol(time.substr(0, 2));
    THiSC[1] = stol(time.substr(3, 2));
    THiSC[2] = stol(time.substr(6, 2));
    THiSC[3] = ns;
    return 0;
}
long int GetDT(Event &evt, int THiSC[])
{
    long int dT0 = (evt.MuStn(0).Traw()[0] - THiSC[0]) * 3600 +(evt.MuStn(0).Traw()[1] - THiSC[1]) * 60 + (evt.MuStn(0).Traw()[2] - THiSC[2]);
    long int dT = dT0 * 1000000000 + evt.MuStn(0).Traw()[3] + evt.MuStn(0).Delay() - THiSC[3];
    return dT;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Write to root-file    events     from Grande coinciding with HiSCORE events
int SelectGrHiscEvts(const char *grDat = GRANDE_DATA_ADDRESS, const char* hiscDat = HiSCORE_DATA_ADDRESS, const char* fileToWrite = OUTPUT_FILE_NAME)
{
    vector<Event> evts, evtsNew;
    ReadTreeOfEvts(grDat, evts);

    long int dT, dTlast;
    string line;
    ifstream hiscore(hiscDat);
    if (! hiscore.is_open()) return -1;
    int ii = 0;
    while(getline(hiscore, line))
    {
        //ii++;
        int k = 0;
        int THiSC[4];
        GetHiSCtime(':', line, k, THiSC);
        dTlast = DT_WINDOW;
        for(int i = 0; i < (int)evts.size(); i++)
        {
            dT = GetDT(evts[i], THiSC);

            if(abs(dT) < DT_WINDOW || (abs(dT) >= DT_WINDOW && abs(dTlast) < DT_WINDOW))
            {
                if((abs(dTlast) < abs(dT) && dTlast != DT_WINDOW) || (abs(dT) >= DT_WINDOW && abs(dTlast) < DT_WINDOW))
                {
                    i--;
                    dT = dTlast;
                    evts[i].SetEvtTimeHiSC(THiSC);
                    evts[i].SetMatchHiSC(1);//matching Gr and HiSC
                    WriteEASdirection(line, k, evts[i]);
                    WriteEAScenter(line, evts[i], k);
                    ii++;

                    //evtsNew.push_back(evts[i]);
                    h0->Fill(dT/100.);
                    i++;
                    dTlast = DT_WINDOW;
                    break;
                }
                else dTlast = dT;
            }
            else if(abs(dT) > DT_WINDOW && dT > 0)
            {
               // evtsNew.push_back(evts[i]);
                break;
            }
        }
    }
    hiscore.close();
    WriteToFile(evts, fileToWrite);
    cout << ii <<" matched evts"<<endl;
    auto c1 = new TCanvas("c1", "");
    h0 -> Draw();
    return 0;
}
