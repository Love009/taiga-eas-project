#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;

const char *FILE_TO_RECON = "analysis_results/Gr+Hisc261119.root";
bool DEBUG = false;
static const int    N_STN_MIN =   3;
static const int    N_STN_TOT =  19;
Event evt;

TH1D *h0   = new TH1D("h0",   "Hist to fit", 20,  0,  20);
TH2D *h1   = new TH2D("h1", "dT_exp-dT_theory [ns], st_i-st_1", 19, 1, 19, 2000, -1000, 1000);
TH1D *h2   = new TH1D("h2", "dT_exp-dT_theory [ns], 2 neighbour stations", 2000, -1000, 1000);
TH1D *dt_hist   = new TH1D("dt_hist",   "dT_exp - dT_theory [ns]", 1000, -500, 500);
TH1D *psi = new TH1D("psi", "psi angle, degrees",   15*10, 0, 15); //angle between reconstructed EAS vector and generated (or reconstructed by HISCORE) one
//TH1D *psi_ = new TH1D("psi_", "1 - cos(psi angle)^2",   70*4, 0, 0.067668);
TH1D *psi_ = new TH1D("psi_", "sin^2(psi angle / 4)",   70*4, 0, 0.004277);
TH1D *dphi = new TH1D("dphi", "delta phi, degrees",   40*2, -20, 20);
TH1D *dthe = new TH1D("dthe", "delta theta, degrees", 40*2, -20, 20);
TH1D *chi2 = new TH1D("chi2", "chi2", 110, -10, 100);
TH1D *nTrigSts = new TH1D("nTrigSts",   "Number of stations triggered in event", (N_STN_TOT ),  1,  N_STN_TOT);
TH1D *thetaDistr = new TH1D("thetaDistr", "theta, degrees", 190,  -5,  180);
TH1D *phiDistr = new TH1D("phiDistr", "phi, degrees", 370,  -5,  365);
TH1D *thetaDistrRec = new TH1D("thetaDistrRec", "theta, degrees", 190,  -5,  180);
TH1D *phiDistrRec = new TH1D("phiDistrRec", "phi, degrees", 370,  -5,  365);

//TGraph *gr = new TGraph();
//int grSize = 0;


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
    Double_t val = t0 + (n1*evt.MuStn(ist).X() + n2*evt.MuStn(ist).Y() + n3*evt.MuStn(ist).Z())/0.3;

    return val;
}

Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t arg = 0;
    Double_t xx = x[0];
    Double_t fitval = par[0] + par[2] * exp(xx * par[1]);
    return fitval;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create and fit histogram
void FindDirection(Event& evnt)
{
    double xl = 0.0;
    double xr = evnt.NStations();//number of points of histogram
    TF1  *f1 = new TF1("Plane", PlaneFit, xl, xr, 3);
    h0 -> Reset();
    double T_mean = 0;
    for (int i=0; i < evnt.NStations(); i++) { // Filling the histogram to fit
        double t = evnt.MuStn(i).T()-evnt.MuStn(0).T();
        T_mean += t;
        double e = evnt.MuStn(i).Terr();
        h0 -> SetBinContent(i+1,t);
        h0 -> SetBinError(i+1,e*10./3.);
    }
    T_mean = T_mean / evnt.NStations();
    double T_min = evnt.MuStn(0).T(), T_max = evnt.MuStn(0).T();
    int stIDmin = 0;
    double x_min = evnt.MuStn(0).X(), x_max = evnt.MuStn(0).X(), y_min = evnt.MuStn(0).Y(), y_max = evnt.MuStn(0).Y();
    for (int i=0; i < evnt.NStations(); i++) {
        if(evnt.MuStn(i).T() < T_min)
        {
            T_min = evnt.MuStn(i).T();
            stIDmin = i;
            x_min =evnt.MuStn(i).X();
            y_min =evnt.MuStn(i).Y();
        }
        if(evnt.MuStn(i).T() > T_max)
        {
            T_max = evnt.MuStn(i).T();
            x_max =evnt.MuStn(i).X();
            y_max =evnt.MuStn(i).Y();
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
    evnt.SetThetaRec(the);
    evnt.SetPhiRec(phi);
    evnt.SetChi2(chi);
    int *TimeCust;
    TimeCust = evnt.MuStn(stIDmin).Trow();
    TimeCust[3] = (int) t0;
    evnt.SetEvtTime(TimeCust);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FindDT(Event &evt, int num1, int num2)
{
    double dt_exp = double(evt.MuStn(num1).T() - evt.MuStn(num2).T());
    double x, y, z;
    x = evt.MuStn(num1).X() - evt.MuStn(num2).X();
    y = evt.MuStn(num1).Y() - evt.MuStn(num2).Y();
    z = evt.MuStn(num1).Z() - evt.MuStn(num2).Z();
    double n1 = cos(evt.PhiRec())*sin(evt.ThetaRec());
    double n2 = sin(evt.PhiRec())*sin(evt.ThetaRec());
    double n3 = cos(evt.ThetaRec());
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
    auto c1 = new TCanvas("c1", "");
    c1->Divide(1,2);
    c1->cd(1);
    //gr -> Draw("AP*");

    psi_ -> Scale(1. / psi_->GetEntries());
    TF1 *func = new TF1("fit", fitf, 0., 0.004, 3);
    func -> SetParameters(0., -1000., 1.);
    func -> FixParameter(2, psi_ ->GetMaximum());
   // func -> SetParameter(2, psi_ ->GetMaximum());
    func -> SetParLimits(0, 0., 1.);
    func -> SetParLimits(1, -1000000, 0.);
   // func -> SetParLimits(2, 0., 1.);
    psi_ -> Fit("fit", "", "", 0., 0.004);
    c1->SetLogy();
    psi_->Draw();
    cout << psi_ ->GetMaximum() << "real\t vs" << func -> GetParameter(0)+ func -> GetParameter(2) << endl;

    gStyle->SetOptStat(111);//output without RMS
    c1->cd(2);
    //auto c1 = new TCanvas("c1", "");
    //c1->SetLogy();
    psi->Draw();
    
   /* auto c2 = new TCanvas("c2", "");
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
    h2->Draw();*/
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FillAnglHists(Event &evt, int &counterOfEvts_chi2)
{
    int toOptimize = 1;
    do{
        FindDirection(evt);
        double d_psi = sin(evt.ThetaRec())*cos(evt.PhiRec())*sin(evt.ThetaGen())*cos(evt.PhiGen())
                       + sin(evt.ThetaRec())*sin(evt.PhiRec())*sin(evt.ThetaGen())*sin(evt.PhiGen())
                       + cos(evt.ThetaRec())*cos(evt.ThetaGen());
        double dAngle = evt.PhiRec() - evt.PhiGen();
        if(dAngle >= TMath::Pi()) dAngle = 2 *TMath::Pi() - dAngle;
        if(dAngle <= -TMath::Pi()) dAngle = 2 * TMath::Pi() + dAngle;
        if(evt.Chi2() < 10)
        {
            toOptimize = 0;
            counterOfEvts_chi2++;
            if (DEBUG)    cout  << "Reconstructed angle Theta = " << evt.ThetaRec() * 57.3 <<  " : " << evt.ThetaGen() * 57.3 << "\n"  << "Reconstructed angle Phi   = "   << evt.PhiRec() * 57.3  <<  " : " <<evt.PhiGen() * 57.3  << "\n"  << "Chi2 of the fit           = "            << evt.Chi2()     <<"\n" << endl;

            //grSize++;
            //gr->Expand(grSize);
            //gr->SetPoint(grSize - 1, TMath::ACos(d_psi), 1 - pow(TMath::Cos(TMath::ACos(d_psi)),2));
            
            psi -> Fill(57.3 * TMath::ACos(d_psi));
            psi_ -> Fill(pow(TMath::Sin(TMath::ACos(d_psi) / 4.),2));
            double angle = evt.ThetaRec();
            if(evt.ThetaRec() < 0) angle = angle * (-1);
            thetaDistr->Fill(57.3*evt.ThetaGen());
            thetaDistrRec->Fill(57.3*angle);

            dthe -> Fill(57.3 * (angle - evt.ThetaGen()));
            dphi -> Fill(57.3 * dAngle*TMath::Sin(evt.ThetaRec()));

            angle = evt.PhiGen();
            angle = CorrectAngleVal(angle);
            phiDistr->Fill(57.3*angle);
            angle = evt.PhiRec();
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
int FillTimeHists(Event &evt)
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
        if(evts[ievt].NStations() > 2) nTrigSts->Fill(evts[ievt].NStations());
        if(evts[ievt].NStations() >= N_STN_MIN)
        {
            counterOfEvts++;
            evt.Reset();
            evt = evts[ievt];
            if(evts[ievt].MatchHiSC() == 1){
                FillAnglHists(evt, counterOfEvts_chi2);
                FillTimeHists(evt);
            }

        }
    }
    cout << "Evts num before chi2 " << counterOfEvts << ", after chi2 " << counterOfEvts_chi2 << endl;
    WriteToFile(evts, FileName);
    if(draw) DrawHists();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Some custom functions for quick analysis:
int main1()
{
    gErrorIgnoreLevel = kError; //without warnings about TCanvas deleting
    Recon("analysis_results/Gr+Hisc130220.root", 0);
    Recon("analysis_results/Gr+Hisc150220.root", 0);
    Recon("analysis_results/Gr+Hisc160220.root", 0);
    Recon("analysis_results/Gr+Hisc200220.root", 0);
    Recon("analysis_results/Gr+Hisc210220.root", 0);
    Recon("analysis_results/Gr+Hisc250220.root", 0);
    Recon("analysis_results/Gr+Hisc270220.root");
    return 0;
}

