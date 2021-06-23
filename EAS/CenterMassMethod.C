#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;

static const int    N_STN_MIN =   4;
static const int    N_STN_TOT =  19;
static const int    INTEGRATION_INTERVAL = 60;
double coord[2];

vector<double> MeanCH0;
vector<double> MeanCH2;
vector<double> SigmaCH0;
vector<double> SigmaCH2;
const char *INITIALIZATION_FILE_NAME = "analysis_results/pedestals/pedestal241119.txt";
const char *FILE_TO_READ = "analysis_results/for_EAS_center/Gr+Hisc241119Recon.root";

TH1D *dx = new TH1D("dx", "delta x, m",   120, -200, 200);
TH1D *dy = new TH1D("dy", "delta y, m",   120, -200, 200);
TH2D *GenRecX  = new TH2D("GenRecX", "Gen x and rec x, m", 1000, -500, 500, 1000, -500, 500);
TH2D *GenRecY  = new TH2D("GenRecY", "Gen y and rec y, m", 1000, -500, 500, 1000, -500, 500);
TH2D *EAScenter = new TH2D("EAScenter", "EAS center, m", 1000., -500., 500., 1000., -500., 500.);

const int TIME_WINDOW = 5; // *5 ns
const int LOW_SIG_CHECK_INTERVAL = 30; // counts of 200 total ADS counts
//===============================================================================================//
//Initualize station channel pedestals
void Initialize(const char *InitFileName = INITIALIZATION_FILE_NAME)
{
    double meanCH0, meanCH2, sigmaCH0, sigmaCH2;
    string  title, rest;
    int i = 1;
    ifstream in(InitFileName);
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Find a signal start using linear fit (the least square method)
int FindSignal(vector<int> ADCsig, int stID, int chanID, double sigma, int threshold, int lowLimit, int window = TIME_WINDOW)
{
    int sigTime = -1., flag = 0;
    double mean = MeanCH0[stID - 1];
    if(chanID == 2) mean = MeanCH2[stID - 1];
    for(int k = 0; k < (int) ADCsig.size(); k++)
    {
        int delta = ADCsig[k] -  mean;
        if(delta <= lowLimit * sigma && k <= LOW_SIG_CHECK_INTERVAL)
        {
            cout << "Signal ampl is lower than " << lowLimit * sigma << "\tStation_number = " << stID << "\t channel " << chanID << "\n\n";
            sigTime = -2.;
            break;
        }
        else if(flag >= window)
        {
            TArrayD x(2*window), y(window);
            TMatrixD A(window, 2), A_t(2, window),Y(window, 1), Pars(2, 1);
            for(int i = 0; i < window; i++)
            {
                y[i] = ADCsig[k - 1 - window + i] - mean;
                x[2*i] = k - 1 - window + i;
                x[2*i+1] = 1;
            }
            A.SetMatrixArray(x.GetArray());
            Y.SetMatrixArray(y.GetArray());
            A_t.Transpose(A);
            TMatrixD M = A_t*A;
            M.Invert();
            Pars = M * A_t * Y;
            double *pars = Pars.GetMatrixArray();
            sigTime = -pars[1] / pars[0];
            break;
        }
        else if(delta >= threshold * sigma) flag++;
        else flag = 0;
    }
    return sigTime;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FindSignalStart(Event &evt, int j, int &sigTime0, int &sigTime2, int threshold = 5, int lowLimit = -10)
{
    int res;
    int stID = evt.MuStn(j).ID();
    sigTime0 = FindSignal(evt.MuStn(j).ADC0(), stID, 0, SigmaCH0[stID - 1], threshold, lowLimit);
    sigTime2 = FindSignal(evt.MuStn(j).ADC2(), stID, 2, SigmaCH2[stID - 1], threshold, lowLimit);
    if(sigTime0 >= 0 && sigTime2 >= 0) res = 1;
    else res = 0;
    return res;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FindIntegrals(Event &evt, double *integralValues)
{
    int res = 1;
    for(int j = 0; j < evt.NStations(); j++)
    {
        double integral0 = 0, integral2 = 0;
        int flag0 = 0, flag2 = 0;
        int max0 = evt.MuStn(j).ADC0()[0];
        int max2 = evt.MuStn(j).ADC2()[0];
        int sigStart0 = 0, sigStart2 = 0;

        if(FindSignalStart(evt, j, sigStart0, sigStart2))
        {
            for (int k = 0; k < (int) evt.MuStn(j).ADC0().size(); k++)
            {
                if (k >= sigStart0 && (k - sigStart0) <= INTEGRATION_INTERVAL) integral0 += evt.MuStn(j).ADC0()[k] -
                                                                                            MeanCH0[evt.MuStn(j).ID() -
                                                                                                    1];
                if (k >= sigStart2 && (k - sigStart2) <= INTEGRATION_INTERVAL) integral2 += evt.MuStn(j).ADC2()[k] -
                                                                                            MeanCH2[evt.MuStn(j).ID() -
                                                                                                    1];

                //Check saturation
                if (evt.MuStn(j).ADC0()[k] > max0) max0 = evt.MuStn(j).ADC0()[k];
                if (evt.MuStn(j).ADC0()[k] == max0 &&
                    evt.MuStn(j).ADC0()[k] > 1.9 * MeanCH0[evt.MuStn(j).ID() - 1])
                    flag0++;
                if (evt.MuStn(j).ADC2()[k] > max2) max2 = evt.MuStn(j).ADC2()[k];
                if (evt.MuStn(j).ADC2()[k] == max2 &&
                    evt.MuStn(j).ADC2()[k] > 1.9 * MeanCH2[evt.MuStn(j).ID() - 1])
                    flag2++;
                if (flag0 == 3 || flag2 == 3)
                {
                    coord[0] = -10000;
                    coord[1] = -10000;
                    res = -1;
                    break;
                }
                //end of saturation check
            }
            integralValues[j] = (integral0 + integral2);
        }
        else res = -1;
    }
    return res;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//15.02.2021 (Integral method)
//Find the max signal and return the station number ---- from 0 to (NStations-1)---- in the event evt (!not the station ID!)
//USING INTEGRALS
int FindCoords(Event evt)
{
    int n = evt.NStations();
    if(n < N_STN_MIN)
    {
        coord[0] = -10000;
        coord[1] = -10000;
        return -1;
    }
    double integralValues[n];
    if(FindIntegrals(evt, integralValues))
    {
        double sum = 0;
        coord[0] = 0;
        coord[1] = 0;
        for (int i = 0; i < n; i++) sum += integralValues[i];
        for (int i = 0; i < n; i++)
        {
            coord[0] += integralValues[i] * evt.MuStn(i).X() / sum;
            coord[1] += integralValues[i] * evt.MuStn(i).Y() / sum;
        }
    }
    else return -1;
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CheckCenter400(Event &evt)
{
    if(sqrt(evt.EAScentersGen()[0] * evt.EAScentersGen()[0] + evt.EAScentersGen()[1] * evt.EAScentersGen()[1]) <= 400.) return 1;
    else return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CheckHiSCAperture(Event &evt)
{
    double HiSC_the = 20. / (180. / TMath::Pi());
    double HiSC_view_size = 30. / (180. / TMath::Pi());

    double the = evt.ThetaRec() , phi = evt.PhiRec();
    if(TMath::ACos(sin(HiSC_the) * sin(the) * cos(phi) + cos(HiSC_the) * cos(the)) <= HiSC_view_size) return 1;
    else return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int WriteToFile(vector<Event> &evts, const char* FileName)
{
    string fileToWrite(FileName);
    fileToWrite = fileToWrite.substr(0, (int)strlen(FileName) - 5);
    fileToWrite += "Center.root";
    TFile f(fileToWrite.c_str(),"recreate");
    TTree Muon("Muon","");
    gROOT->ProcessLine("#include <vector>");
    Event *event = new Event();
    Muon.Branch("events", &event);
    Muon.SetAutoSave(1000000);
    for(int i = 0; i < (int)evts.size(); i++)
    {
        event = & evts[i];
        Muon.Fill();
    }
    f.Write();
    f.Close();
    cout << "I have written it!" << endl;
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Draw_()
{
    auto c1 = new TCanvas("c1", "");
    gStyle->SetOptStat(1111);
    c1->Divide(2, 1);
    c1->cd(1);
    dx->GetXaxis()->SetTitle("dx, m");
    dx -> Fit("gaus", "", "", -110,110);
    dx->Draw();
    c1->cd(2);
    dy->GetXaxis()->SetTitle("dy, m");
    dy -> Fit("gaus", "", "", -110,110);
    dy->Draw();

    auto c2 = new TCanvas("c2", "");
    c2->Divide(2, 1);
    c2->cd(1);
    GenRecX->SetMarkerStyle(7);
    GenRecX->GetXaxis()->SetTitle("x HiSCORE, m");
    GenRecX->GetYaxis()->SetTitle("x RecGrande, m");
    GenRecX->GetYaxis()->SetTitleOffset(1.5);
    GenRecX->Draw();
    c2->cd(2);
    GenRecY->SetMarkerStyle(7);
    GenRecY->GetXaxis()->SetTitle("y HiSCORE, m");
    GenRecY->GetYaxis()->SetTitle("y RecGrande, m");
    GenRecY->GetYaxis()->SetTitleOffset(1.5);
    GenRecY->Draw();

    auto c17 = new TCanvas("c17", "");
    //EAScenter->SetMarkerStyle(7);
    // EAScenter->SetMarkerSize(1);
    EAScenter->Draw();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FindEASCenter(const char *FileName = FILE_TO_READ, int draw = 1, const char *InitFileName = INITIALIZATION_FILE_NAME)
{
    Initialize(InitFileName);
    vector<Event> evts;
    ReadTreeOfEvts(FileName, evts);
    int counter = 0;
    for(int i = 0; i < (int)evts.size(); i++)
    {
        int flag400 = CheckCenter400(evts[i]);
        int flagHiSCaperture = CheckHiSCAperture(evts[i]);
        if(evts[i].MatchHiSC() == 1 && evts[i].Chi2() < 10 && flag400 == 1 && flagHiSCaperture == 1 && evts[i].NStations() == N_STN_MIN) 
        {
            int res = FindCoords(evts[i]);
            if (res == -1)
            {
               coord[0] = -10000;
               coord[1] = -10000;
               evts[i].SetEAScenterRec(coord);
               counter++;
            }
            else
            {
                double *xyGen;
                xyGen = evts[i].EAScentersGen();
                //cout << "Grande Hisc: X = " << coord[0]<< "\t" << xyGen[0]<< "\tY = " << coord[1] << "\t" << xyGen[1] << endl;
                dx->Fill(coord[0] - xyGen[0]);
                dy->Fill(coord[1] + xyGen[1]);
                GenRecX->Fill(xyGen[0], coord[0]);
                GenRecY->Fill(-xyGen[1], coord[1]);
                EAScenter -> Fill(coord[0], coord[1]);
                //EAScenter -> Fill(xyGen[0], -xyGen[1]);
                evts[i].SetEAScenterRec(coord);
            }
        }
    }
    WriteToFile(evts, FileName);
    if(draw) Draw_();
    return 0;
}

int main3()
{
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc241119Recon.root", 0, "analysis_results/pedestals/pedestal241119.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc241219Recon.root", 0, "analysis_results/pedestals/pedestal241219.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc291219Recon.root", 0, "analysis_results/pedestals/pedestal291219.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc020120Recon.root", 0, "analysis_results/pedestals/pedestal020120.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc221119Recon.root", 0, "analysis_results/pedestals/pedestal221119.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc261119Recon.root", 0, "analysis_results/pedestals/pedestal261119.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc271119Recon.root", 0, "analysis_results/pedestals/pedestal271119.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc301219Recon.root", 0, "analysis_results/pedestals/pedestal301219.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc150220Recon.root", 0, "analysis_results/pedestals/pedestal150220.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc160220Recon.root", 0, "analysis_results/pedestals/pedestal160220.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc200220Recon.root", 0, "analysis_results/pedestals/pedestal200220.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc250220Recon.root", 0, "analysis_results/pedestals/pedestal250220.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc270220Recon.root", 0, "analysis_results/pedestals/pedestal270220.txt");
    FindEASCenter("analysis_results/for_EAS_center/Gr+Hisc210220Recon.root", 1, "analysis_results/pedestals/pedestal210220.txt");
    return 0;
}
