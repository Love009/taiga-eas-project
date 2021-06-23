#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;

const char *FILE_TO_RECON = "analysis_results/Gr+Hisc261119.root";
bool DEBUG = false;
static const int    N_STN_MIN =   3;
static const int    N_STN_TOT =  19;
EventAnalysis evt;

TH1D *h0   = new TH1D("h0",   "Hist to fit", 20,  0,  20);
TH2D *h1   = new TH2D("h1", "dT_exp-dT_theory [ns], st_i-st_1", 19, 1, 19, 2000, -1000, 1000);
TH1D *h2   = new TH1D("h2", "dT_exp-dT_theory [ns], 2 neighbour stations", 2000, -1000, 1000);
TH1D *dt_hist   = new TH1D("dt_hist",   "dT_exp - dT_theory [ns]", 1000, -500, 500);
TH1D *psi = new TH1D("psi", "psi angle, degrees",   34*2, -2, 15); //angle between reconstructed EAS vector and generated (or reconstructed by HISCORE) one
TH1D *dphi = new TH1D("dphi", "delta phi, degrees",   40*2, -20, 20);
TH1D *dthe = new TH1D("dthe", "delta theta, degrees", 40*2, -20, 20);
TH1D *chi2 = new TH1D("chi2", "chi2", 110, -10, 100);
TH1D *nTrigSts = new TH1D("nTrigSts",   "Number of stations triggered in event", (N_STN_TOT ),  1,  N_STN_TOT);
TH1D *thetaDistr = new TH1D("thetaDistr", "theta, degrees", 190,  -5,  180);
TH1D *phiDistr = new TH1D("phiDistr", "phi, degrees", 370,  -5,  365);
TH1D *thetaDistrRec = new TH1D("thetaDistrRec", "theta, degrees", 190,  -5,  180);
TH1D *phiDistrRec = new TH1D("phiDistrRec", "phi, degrees", 370,  -5,  365);
//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Fit function
Double_t PlaneFit(Double_t *x, Double_t *par)
{
    Float_t xx = x[0];
    int ist = int(xx);

    Double_t the = par[0];
    Double_t phi = par[1];
    Double_t t0  = par[2];

    Double_t n1 = sin(the)*cos(phi);
    Double_t n2 = sin(the)*sin(phi);
    Double_t n3 = cos(the);
    Double_t val = t0 + (n1*evt.MuStn(ist).X()[0] + n2*evt.MuStn(ist).Y()[0] + n3*evt.MuStn(ist).Z()[0])/0.3;

    return val;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create and fit histogram
void FindDirection(EventAnalysis& evnt)
{
    double xl = 0.0;
    double xr = evnt.NStations();//number of points of histogram
    TF1  *f1 = new TF1("Plane", PlaneFit, xl, xr, 3);
    h0 -> Reset();
    double T_mean = 0;
    for (int i=0; i < evnt.NStations(); i++) { // Filling the histogram to fit
        double t = double(evnt.MuStn(i).Traw()[3] + evnt.MuStn(i).Delay() - evnt.MuStn(0).Traw()[3] - evnt.MuStn(0).Delay());
        T_mean += t;
        double e = 15;//evnt.MuStn(i).Terr();
        h0 -> SetBinContent(i+1,t);
        h0 -> SetBinError(i+1,e);
    }
    T_mean = T_mean / evnt.NStations();
    double T_min = double(evnt.MuStn(0).Traw()[3] + evnt.MuStn(0).Delay()), T_max = double(evnt.MuStn(0).Traw()[3] + evnt.MuStn(0).Delay());
    int stIDmin = 0;
    double x_min = evnt.MuStn(0).X()[0], x_max = evnt.MuStn(0).X()[0], y_min = evnt.MuStn(0).Y()[0], y_max = evnt.MuStn(0).Y()[0];
    for (int i=0; i < evnt.NStations(); i++) {
        if(double(evnt.MuStn(i).Traw()[3] + evnt.MuStn(i).Delay()) < T_min)
        {
            T_min = double(evnt.MuStn(i).Traw()[3] + evnt.MuStn(i).Delay());
            stIDmin = i;
            x_min =evnt.MuStn(i).X()[0];
            y_min =evnt.MuStn(i).Y()[0];
        }
        if(double(evnt.MuStn(i).Traw()[3] + evnt.MuStn(i).Delay()) > T_max)
        {
            T_max = double(evnt.MuStn(i).Traw()[3] + evnt.MuStn(i).Delay());
            x_max =evnt.MuStn(i).X()[0];
            y_max =evnt.MuStn(i).Y()[0];
        }
    }
    if(DEBUG)cout << "phi_0 = " << 57.3 * TMath::ATan((y_max-y_min)/(x_max-x_min)) << "(x,y)_min " << x_min << "  " << y_min<< "(x,y)_max " << x_max << " " <<y_max << " T_min, T_max = " << T_min << "   " << T_max <<endl;
    f1->SetParameters(20./57.3, TMath::ATan((y_max-y_min)/(x_max-x_min)), T_mean);
    f1->SetParLimits(0, -1.58, 1.58);
    f1->SetParLimits(1, -6.30, 6.30); //TMath::Pi() * 2);
    f1->SetParLimits(2,  -20000, 20000);
    f1->SetLineStyle(3);

//Fitting the histogram
    evt = evnt;
    h0->Fit("Plane", "Q0", "", xl, xr); // "Q" - no printing, "P" - Pearson chi-square method
    double the = f1 -> GetParameter(0);
    double phi = f1 -> GetParameter(1);
    double t0  = f1 -> GetParameter(2);
    double chi = f1 -> GetChisquare(); // chi = hi_squared
    double prp = f1 -> GetProb();
    int    ndf = f1 -> GetNDF();
//Writing results
    evnt.SetThetaGr(the);
    evnt.SetPhiGr(phi);
    evnt.SetChi2(chi);
    /*long int *TimeCust;
    TimeCust = evnt.MuStn(stIDmin).Traw();
    TimeCust[3] = (long int) t0;
    evnt.SetEvtTime(TimeCust);*/
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FindDT(EventAnalysis &evt, int num1, int num2)
{
    double temp1, temp2, temp3, temp;
    temp1 = double((evt.MuStn(num1).Traw()[0] - evt.MuStn(num2).Traw()[0]));
    temp2 = double((evt.MuStn(num1).Traw()[1] - evt.MuStn(num2).Traw()[1]));
    temp3 = double((evt.MuStn(num1).Traw()[2] - evt.MuStn(num2).Traw()[2]));
    temp = temp1 * 3600 + temp2 * 60 + temp3;
    double dt_exp = temp * 1000000000 + double(evt.MuStn(num1).Traw()[3] - evt.MuStn(num2).Traw()[3] + evt.MuStn(num1).Delay() - evt.MuStn(num2).Delay());
    double x, y, z;
    x = evt.MuStn(num1).X()[0] - evt.MuStn(num2).X()[0];
    y = evt.MuStn(num1).Y()[0] - evt.MuStn(num2).Y()[0];
    z = evt.MuStn(num1).Z()[0] - evt.MuStn(num2).Z()[0];
    double n1 = cos(evt.PhiGr())*sin(evt.ThetaGr());
    double n2 = sin(evt.PhiGr())*sin(evt.ThetaGr());
    double n3 = cos(evt.ThetaGr());
    double dt_theory = (n1*x+n2*y+n3*z)/0.3;
    double dT = dt_exp - dt_theory;
    return dT;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CorrectAngleVal(double angle)
{
    if(angle >= 2*TMath::Pi()) angle = -2 *TMath::Pi()+angle;
    if(angle <= -2*TMath::Pi()) angle =-(angle + 2 * TMath::Pi() );
    if(angle < 0) angle = 2 * TMath::Pi() + angle;
    return angle;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int WriteToFile(vector<Event> &evts, const char* FileName)
{
    string fileToWrite(FileName);
    fileToWrite = fileToWrite.substr(0, (int)strlen(FileName) - 5);
    fileToWrite += "Recon.root";
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
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int DrawHists()
{
    gStyle->SetOptStat(111);//output without RMS
    auto c1 = new TCanvas("c1", "");
    //c1->SetLogy();
    psi->Draw();
    auto c2 = new TCanvas("c2", "");
    c2->Divide(2, 1);
    c2->cd(1);
    //gPad->SetLogy();
    dphi->Draw();
    c2->cd(2);
    dthe->Draw();
    auto c3 = new TCanvas("c3", "");
    gPad->SetLogy();
    chi2->Draw();
    auto c4 = new TCanvas("c4", "");
    gPad->SetLogy();
    dt_hist->Draw();
    auto c5 = new TCanvas("c5", "");
    h1->SetMarkerStyle(kFullTriangleDown);
    h1->Draw();
    auto c6 = new TCanvas("c6", "");
    gPad->SetLogy();
    h2->Draw();
    auto c10 = new TCanvas("c10", "");
    nTrigSts -> Draw();
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FillAnglHists(EventAnalysis &evt, int &counterOfEvts_chi2)
{
    int toOptimize = 1;
    do{
        FindDirection(evt);
        double d_psi = sin(evt.ThetaGr())*cos(evt.PhiGr())*sin(evt.ThetaHiSC())*cos(evt.PhiHiSC())
                       + sin(evt.ThetaGr())*sin(evt.PhiGr())*sin(evt.ThetaHiSC())*sin(evt.PhiHiSC())
                       + cos(evt.ThetaGr())*cos(evt.ThetaHiSC());
        double dAngle = evt.PhiGr() - evt.PhiHiSC();
        if(dAngle >= TMath::Pi()) dAngle = 2 *TMath::Pi() - dAngle;
        if(dAngle <= -TMath::Pi()) dAngle = 2 * TMath::Pi() + dAngle;
        if(evt.Chi2() < 10)
        {
            toOptimize = 0;
            counterOfEvts_chi2++;
            if (DEBUG)    cout  << "Reconstructed angle Theta = " << evt.ThetaGr() * 57.3 <<  " : " << evt.ThetaHiSC() * 57.3 << "\n"  << "Reconstructed angle Phi   = "   << evt.PhiGr() * 57.3  <<  " : " <<evt.PhiHiSC() * 57.3  << "\n"  << "Chi2 of the fit           = "            << evt.Chi2()     <<"\n" << endl;
            psi -> Fill(57.3 * TMath::ACos(d_psi));
            double angle = evt.ThetaGr();
            if(evt.ThetaGr() < 0) angle = angle * (-1);
            thetaDistr->Fill(57.3*evt.ThetaHiSC());
            thetaDistrRec->Fill(57.3*angle);

            dthe -> Fill(57.3 * (angle - evt.ThetaHiSC()));
            dphi -> Fill(57.3 * dAngle*TMath::Sin(evt.ThetaGr()));

            angle = evt.PhiHiSC();
            angle = CorrectAngleVal(angle);
            phiDistr->Fill(57.3*angle);
            angle = evt.PhiGr();
            angle = CorrectAngleVal(angle);
            phiDistrRec->Fill(57.3*angle);

        }
        else if(evt.Chi2() >= 10 && evt.NStations() > 3)
        {
            double dT, dTlast = 0;
            int badStnNum = 0;
            for(int l = 1; l < evt.NStations(); l++)
            {
                dT = FindDT(evt, l, badStnNum);
                if(dT > dTlast) {badStnNum = l; dTlast = dT;}
            }
            if(dTlast != 0){ evt.DeleteStation(badStnNum); toOptimize = 1; cout << "dT_BAD = (dt_exp - dt_theory) = " << dTlast << endl; }
            else toOptimize = 0;
        }
        else toOptimize = 0;
    }while(toOptimize);
    if(evt.NStations() > 3) chi2 -> Fill(evt.Chi2());
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FillTimeHists(EventAnalysis &evt)
{
    for(int j = 1; j < evt.NStations(); j++)
    {
        int min_j = 0;
        double dT = FindDT(evt, j, min_j);
        if(evt.NStations() >= 3) dt_hist->Fill(dT);
    }
    for(int j = 0; j < evt.NStations(); j++)
    {
        for(int k = 0; k < evt.NStations(); k++)
        {
            if(k != j && evt.NStations() >= 3)
            {
                double dT = FindDT(evt, k, j);
                dt_hist->Fill(dT);
                if(evt.MuStn(j).ID() == 1) h1 -> Fill(evt.MuStn(k).ID(), dT);
                if(-evt.MuStn(j).ID() + evt.MuStn(k).ID() == 1) h2 -> Fill(dT);
            }
        }
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Reconstruct phi and theta angles of EASes and compare with the angles obtainedby HiSCORE (or with generated in ToyMC angles)
void Recon(const char *FileName = FILE_TO_RECON, int draw = 1)
{
    vector<Event> evts;
    ReadTreeOfEvts(FileName, evts);
    int counterOfEvts = 0, counterOfEvts_chi2 = 0;
    for(int ievt = 0; ievt < (int) evts.size(); ievt++)
    {
        /*evt.Reset();
        EventAnalysis evtTemp(evts[ievt]);
        evt = evtTemp;
        cout << evt.PhiHiSC() << "    " << evts[ievt].PhiHiSC() << "    " << evts[ievt].MuStn(0).Traw()[0] <<  "    " << evt.MuStn(0).Traw()[0] << endl;*/
        if(evts[ievt].NStations() > 2 && evts[ievt].MatchHiSC() == 0) nTrigSts->Fill(evts[ievt].NStations());
        if(evts[ievt].NStations() >= N_STN_MIN && evts[ievt].MatchHiSC() == 1)
        {
            counterOfEvts++;
            evt.Reset();
            EventAnalysis evtTemp(evts[ievt]);
            evt = evtTemp;
            FillAnglHists(evt, counterOfEvts_chi2);
            FillTimeHists(evt);
        }
    }
    cout << "Evts num before chi2 " << counterOfEvts << ", after chi2 " << counterOfEvts_chi2 << endl;
   // WriteToFile(evts, FileName);
    if(draw) DrawHists();
}

