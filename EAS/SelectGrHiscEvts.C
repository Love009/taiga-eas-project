#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;

static const int    N_STN_TOT =  19;
static const int    N_STN_MIN =   2;//for systematic errors analysis N_STN_MIN = 2 is enough
const char*         HiSCORE_DATA_ADDRESS   = "HiSCORE_data/2019-2020/out_2611.dat";
const char*         GRANDE_DATA_ADDRESS    = "analysis_results/Gr261119.root";
const char*         OUTPUT_FILE_NAME       = "analysis_results/test.root";

long int    DT_WINDOW = 10000;// ns
long int    DT_WINDOW0 = 1000000000; // ns
TH1D *h00   = new TH1D("h00", "dT_Gr_HiSC", 2000000,  -1000000,  1000000);
int HiSCstID, HiSCstNum;
TH2D *EAScenterHist = new TH2D("EAScenterHist", "EAS center, m", 1000., -500., 500., 1000., -500., 500.);

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
    if (line[k - 2] != ' ' && line[k - 3] != ' ' && line[k - 4] != ' ') {
        centerX = stod(line.substr(k - 4, 4));

    }
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
    evt.SetEAScenterGen(coordArr);

    EAScenterHist->Fill(centerX, -centerY);
    return 0;
}
int WriteEASdirection(string &line, int k, Event &evt)
{
    evt.SetThetaGen(stod(line.substr(k + 20, 5)) / 57.3);
    double pphi = stod(line.substr(k + 26, 6));
    evt.SetPhiGen(2 * TMath::Pi() - pphi / 57.3);//!!!!!!!!!!!!!!!!!!!!!!!!!
    return 0;
}
int GetHiSCtime(const char sep, string line, int &k, int *THiSC)
{
    while(char(line[k]) != sep)  k++;
    int ns = stoi(line.substr(k + 7, 3)) * 1000000 + stoi(line.substr(k + 11, 3)) * 1000 + stoi(line.substr(k + 15, 3));
    string time = line.substr(k - 2, 9) + to_string(ns);
    THiSC[0] = stoi(time.substr(0, 2));
    THiSC[1] = stoi(time.substr(3, 2));
    THiSC[2] = stoi(time.substr(6, 2));
    THiSC[3] = ns;
    HiSCstID = stoi(line.substr(k - 6, 2));
    HiSCstNum = stoi(line.substr(k + 82, 3));
    return 0;
}
long int GetDT(Event &evt, int THiSC[])
{
    long int dT0 = (evt.EvtTime()[0] - THiSC[0]) * 3600 +(evt.EvtTime()[1] - THiSC[1]) * 60 + (evt.EvtTime()[2] - THiSC[2]);
    long int dT = dT0 * 1000000000 + evt.EvtTime()[3] - THiSC[3];
    return dT;
}
int FindTimeShift(const char* hiscDat, vector<Event> &evts)
{
    long int dT, dTlast;
    string line;
    ifstream hiscore(hiscDat);
    if (! hiscore.is_open()) { cout << "no file" << endl; return -1;}
    int ii = 0;
    while(ii != 10000)
    {
        getline(hiscore, line);
        ii++;
        int k = 0;
        int THiSC[4];
        GetHiSCtime(':', line, k, THiSC);
        dTlast = DT_WINDOW0;
        for(int i = 0; i < (int)evts.size(); i++)
        {
            dT = GetDT(evts[i], THiSC);
            if(abs(dT) < DT_WINDOW0 || (abs(dT) >= DT_WINDOW0 && abs(dTlast) < DT_WINDOW0))
            {
                if((abs(dTlast) < abs(dT) && dTlast != DT_WINDOW0) || (abs(dT) >= DT_WINDOW0 && abs(dTlast) < DT_WINDOW0))
                {
                    dT = dTlast;
                    h00->Fill(dT/1000.);
                    dTlast = DT_WINDOW0;                    
                    break;
                }
                else dTlast = dT;
            }
            else if(abs(dT) > DT_WINDOW0 && dT > 0) break;
        }
    }
    int binmax = h00->GetMaximumBin();
    int shift = h00->GetXaxis()->GetBinCenter(binmax);
    cout << "shift = " << shift*1000 << " ns" << endl;
    hiscore.close();
    return -shift * 1000;
}

int FindMatchedEvts(const char* hiscDat, vector<Event> &evts, int HiSC_CALIBR)
{
    long int dT, dTlast;
    string line;
    ifstream hiscore(hiscDat);
    if (! hiscore.is_open()) return -1;
    int ii = 0;
    TH1D *h0   = new TH1D("h0", "dT_Gr_HiSC, ms", 200,  (-10000 - HiSC_CALIBR)/1000.,  (10000 - HiSC_CALIBR)/1000.);
    int counter600 = 0;
    while(getline(hiscore, line))
    {
        ii++;
        int k = 0;
        int THiSC[4];
        GetHiSCtime(':', line, k, THiSC);
        dTlast = DT_WINDOW;
        for(int i = 0; i < (int)evts.size(); i++)
        {
            dT = GetDT(evts[i], THiSC);
            dT += HiSC_CALIBR;
            if(abs(dT) < DT_WINDOW || (abs(dT) >= DT_WINDOW && abs(dTlast) < DT_WINDOW))
            {
                if((abs(dTlast) < abs(dT) && dTlast != DT_WINDOW) || (abs(dT) >= DT_WINDOW && abs(dTlast) < DT_WINDOW) )
                {
                        i--;
                        dT = dTlast;
                        WriteEASdirection(line, k, evts[i]);
                        WriteEAScenter(line, evts[i], k);
                        evts[i].SetEvtTimeHiSC(THiSC);
                        counter600++;
                        evts[i].SetMatchHiSC(1);
                        dTlast = DT_WINDOW;
                        break;
                }
                else dTlast = dT;
            }
            else if(abs(dT) > DT_WINDOW && dT > 0) break;
        }
    }
    hiscore.close();
    auto c1 = new TCanvas("c1", "");
    h0 -> Draw();
    cout << "counter600 = " << counter600 << endl;
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Write to root-file events     from Grande coinciding with HiSCORE events
int SelectGrHiscEvts(const char *grDat = GRANDE_DATA_ADDRESS, const char* hiscDat = HiSCORE_DATA_ADDRESS, const char* fileToWrite = OUTPUT_FILE_NAME)
{
    double HiSC_CALIBR = 0.;//ns
    vector<Event> evts;
    ReadTreeOfEvts(grDat, evts);
    HiSC_CALIBR = FindTimeShift(hiscDat, evts);
    FindMatchedEvts(hiscDat, evts,HiSC_CALIBR);
    WriteToFile(evts, fileToWrite);
    EAScenterHist-> Draw();
    return 0;
}
