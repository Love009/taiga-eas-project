#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;
#include "TMath.h"

const char *FILE_TO_RECON = "analysis_results/Gr+Hisc261119.root";
bool DEBUG = false;
static const int    N_STN_MIN =  3;
static const int    N_STN_TOT =  19;
static const int    errTemp =  10.;//ns
double dt_syst[N_STN_TOT];


Event evt;

TH1D *h0   = new TH1D("h0",   "Hist to fit", 20,  0,  20);
TH1D *psi = new TH1D("psi", "psi angle, degrees",   100., 0., 10.);
TH1D *cosPsi = new TH1D("cosPsi", "1-cos(psi)", 150., 0., 0.0152);
int n_bins = 25. * 2.5;
TH2D *theta_phi_rec = new TH2D("theta_phi_rec", "theta phi rec, degrees", n_bins,  -3.,  3., n_bins,  -3.,  3.);
TH2D *theta_phi_gen = new TH2D("theta_phi_gen", "theta phi gen, degrees", n_bins,  -3.,  3., n_bins,  -3., 3.);
auto c52 = new TCanvas("c52", "");
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create and fit histogram
void FindDirection(Event& evnt)
{
    double xl = 0.0;
    double xr = evnt.NStations();
    TF1  *f1 = new TF1("Plane", PlaneFit, xl, xr, 3);
    h0 -> Reset();
    double T_mean = 0;

    for (int i = 0; i < evnt.NStations(); i++) { // Filling the histogram to fit
        double t = evnt.MuStn(i).T();
        T_mean += t;
        double e = errTemp;
        h0 -> SetBinContent(i+1,t);
        h0 -> SetBinError(i+1,e);
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
    f1->SetParameters(20./180. * TMath::Pi(), TMath::ATan((y_max-y_min)/(x_max-x_min)), T_mean);
   // cout << TMath::ATan((y_max-y_min)/(x_max-x_min)) <<"\t" <<  T_mean << endl;
    f1->SetParLimits(0, -TMath::Pi() / 2., TMath::Pi() / 2.);
    f1->SetParLimits(1, -TMath::Pi() * 2, TMath::Pi() * 2);
    f1->SetParLimits(2,  -10000., 10000);
    f1->SetLineStyle(3);

//Fitting the histogram
    evt = evnt;
    h0->Fit("Plane", "Q0", "", xl, xr); // "Q" - no printing, "P" - Pearson chi-square method
    double the = f1 -> GetParameter(0);
    double phi = f1 -> GetParameter(1);
    if(phi < 0. && the > 0.) phi += 2 * TMath::Pi();
    if(the < 0.) {the = -the; phi += TMath::Pi();}
    if(phi < 0.) phi += 2 * TMath::Pi();
    if(phi > TMath::Pi() * 2) phi -= 2 * TMath::Pi();
    double t0  = f1 -> GetParameter(2);
    double chi = f1 -> GetChisquare(); // chi = hi_squared
    double prp = f1 -> GetProb();
    int    ndf = f1 -> GetNDF();
//Writing results
    evnt.SetThetaRec(the);
    evnt.SetPhiRec(phi);
    evnt.SetChi2(chi);
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
double CalcMean(TH2D *hist, int n)
{
    double mean = 0;
    int counter = 0;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            double tmp0 = hist -> GetBinContent(i, j);
            if(tmp0 != 0){ mean += tmp0; counter++;}
        }
    }
    mean = mean / counter;
    return mean;
}
//-------------------------------------------------------------------------------------
double CalcRMS(TH2D *hist, double mean, int n)
{
    double rms = 0;
    int counter = 0;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            double tmp0 = hist -> GetBinContent(i, j);
            if(tmp0 != 0) {rms += (tmp0 - mean) * (tmp0 - mean); counter++;}
        }
    }
    rms = sqrt(rms / counter);
    return rms;
}
//-------------------------------------------------------------------------------------
int ShiftHistContent(TH2D *hist, TH2D *hist2, double mean)
{
    int n = hist->GetXaxis()->GetNbins();

    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            int bin = hist -> GetBin(i, j);
            double tmp = (hist -> GetBinContent(bin)) - mean;
            hist2 -> SetBinContent(bin, tmp);
        }
    }
    return 0;
}
//-------------------------------------------------------------------------------------
int CheckBinning(TH2D *hist, double mean, int n)
{
    int flag = 1;
    double more = 0, less = 0;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            int bin = hist -> GetBin(i, j);
            double tmp = (hist -> GetBinContent(bin)) - mean;
            if(tmp < 0) less++;
            else more++;
        }
    }
    cout << "l/m = " << less / more << endl;
    if(less / more >= 0.8 && less / more <= 1.2) flag = 0;
    cout << "less = " << less << "  more = " << more << endl;
    return flag;
}
//-------------------------------------------------------------------------------------
int ScaleToSigma(TH2D *theta_phi_gen, TH2D *theta_phi_rec)
{
    int flag1 = 0, flag2 = 0;
    int rebin_num = 0;
    double mean_rec, mean_gen, rms_rec, rms_gen;
    int n;
    do {
        if(flag1 == 1 || flag2 == 1) {
            theta_phi_gen->Rebin2D(2,2);
            theta_phi_rec->Rebin2D(2,2);
            rebin_num++;
            flag1 = 0;
            flag2 = 0;
        }
        n = theta_phi_gen->GetXaxis()->GetNbins();
        cout <<"rebin = "<< rebin_num << "  n=" <<n <<endl;
        mean_gen = CalcMean(theta_phi_gen, n);
        mean_rec = CalcMean(theta_phi_rec, n);
        rms_gen = CalcRMS(theta_phi_gen, mean_gen, n);
        rms_rec = CalcRMS(theta_phi_rec, mean_rec, n);
        //if(rebin_num <= 1) flag1 = CheckBinning(theta_phi_gen, mean_gen, n);
        //if(rebin_num <= 1) flag2 = CheckBinning(theta_phi_rec, mean_rec, n);
        cout << "mean_gen = " << mean_gen << "  rms_gen = " << rms_gen << endl;
        cout << "mean_rec = " << mean_rec << "  rms_rec = " << rms_rec << "\n" <<endl;
    }while(flag1 == 1 || flag2 == 1);

    n = theta_phi_gen->GetXaxis()->GetNbins();
    TH2D *theta_phi_rec2 = new TH2D("theta_phi_rec2", "theta_phi_rec2", n, -6., 6., n, -6., 6.);
    TH2D *theta_phi_gen2 = new TH2D("theta_phi_gen2", "theta_phi_gen2", n, -6., 6., n, -6., 6.);
    auto c53 = new TCanvas("c53", "");
    c53->Divide(2,1);
    ShiftHistContent(theta_phi_gen, theta_phi_gen2, mean_gen);
    ShiftHistContent(theta_phi_rec, theta_phi_rec2, mean_rec);
    Double_t scale_rec = 1. / rms_rec;
    Double_t scale_gen = 1. / rms_gen;
    theta_phi_rec2->Scale(scale_rec);
    theta_phi_gen2->Scale(scale_gen);
    c53->cd(1);
    theta_phi_rec2 -> GetXaxis() -> SetTitle("(PhiRec - PhiMoon)*Sin(ThetaMoon)");
    theta_phi_rec2 -> GetYaxis() -> SetTitle("ThetaRec - ThetaMoon");
    theta_phi_rec2 -> Draw("CONT4Z");
    c53->cd(2);
    theta_phi_gen2 -> GetXaxis() -> SetTitle("(PhiGen - PhiMoon)*Sin(ThetaMoon)");
    theta_phi_gen2 -> GetYaxis() -> SetTitle("ThetaGen - ThetaMoon");
    theta_phi_gen2 -> Draw("CONT4Z");
    return 0;
}
//-------------------------------------------------------------------------------------
int DrawHists()
{
    c52->cd(1);
    theta_phi_rec -> GetXaxis() -> SetTitle("(PhiRec - PhiMoon)*Sin(ThetaMoon)");
    theta_phi_rec -> GetYaxis() -> SetTitle("ThetaRec - ThetaMoon");
    theta_phi_rec -> Draw("CONT4Z");
    c52->cd(2);
    theta_phi_gen -> GetXaxis() -> SetTitle("(PhiGen - PhiMoon)*Sin(ThetaMoon)");
    theta_phi_gen -> GetYaxis() -> SetTitle("ThetaGen - ThetaMoon");
    theta_phi_gen -> Draw("CONT4Z");
    auto c3 = new TCanvas("c3", "");
   // c3 -> Divide(1, 2);
   // c3 -> cd(1);
    cosPsi -> SetMinimum(0);
    //cosPsi -> Scale(1. / cosPsi -> GetEntries());

    //cosPsi fitting
    double FitInterval = 0.01;
    TF1 *H1 = new TF1("H1", "pol0", 0., FitInterval);
    H1 -> SetParameter(0, cosPsi -> GetMaximum());
    H1 -> SetParLimits(0, 0., 1. * cosPsi -> GetEntries());
    cosPsi->Fit("H1","L","", 0., FitInterval);
    double chi2H1 = H1 -> GetChisquare();

    TF1 *H2 = new TF1("H2", "[0] *(1 - [1] * exp(- pow(TMath::ACos(1 - x) / (1 - [2]),2)))", 0., FitInterval);
    H2 -> SetParameters(cosPsi -> GetMaximum(), 1. - cosPsi -> GetMinimum() / cosPsi -> GetMaximum(), 5. / 57.3);
    H2 -> SetParLimits(0, 0., 1. * cosPsi -> GetEntries());
    H2 -> SetParLimits(1, 0., 1.);
    H2 -> SetParLimits(2, 0.00001, 180. / 57.3);
    cosPsi -> Fit("H2","L","", 0., FitInterval);
    double chi2H2 = H2 -> GetChisquare();
    cosPsi -> SetMarkerColor(kBlack);
    cosPsi -> SetLineColor(kBlack);
    cosPsi -> Draw("");
    H1 -> SetLineColor(kBlue);
    H1 -> Draw("same");
    H2 -> Draw("same");

    double reliability = sqrt(abs(chi2H1 * chi2H1 - chi2H2 * chi2H2));
    cout <<"chi2H1 = " << chi2H1 << "\tchi2H2 = " << chi2H2 << "\nRel = " << reliability << endl;

    // c3 -> cd(2);
    //psi -> Draw();
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FillAnglHists(Event &evt, int &counterOfEvts_chi2)
{
    int toOptimize = 1;
    do{
        FindDirection(evt);

        double moon_phi = 100. / 57.3;
        double moon_the = 35. / 57.3;

         double d_psi = sin(evt.ThetaRec())*cos(evt.PhiRec())*sin(moon_the)*cos(moon_phi)
                         + sin(evt.ThetaRec())*sin(evt.PhiRec())*sin(moon_the)*sin(moon_phi)
                       + cos(evt.ThetaRec())*cos(moon_the);

        if(evt.Chi2() < 10)
        {
            double temp0 = TMath::ACos(d_psi);
            if(temp0 > TMath::Pi() / 2.) temp0 = TMath::Pi() - temp0;
            cosPsi->Fill(1 - cos(temp0));
            psi -> Fill(180. / TMath::Pi() * temp0);
            toOptimize = 0;
            counterOfEvts_chi2++;
            theta_phi_rec -> Fill(57.3 * (evt.PhiRec() - moon_phi) * TMath::Sin(moon_the), 57.3 * (evt.ThetaRec() - moon_the));
            theta_phi_gen -> Fill(57.3 * (evt.PhiGen() - moon_phi) * TMath::Sin(moon_the), 57.3 * (evt.ThetaGen() - moon_the));
        }
        else if(evt.Chi2() >= 10 && evt.NStations() > 19)
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
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Reconstruct phi and theta angles of EASes and compare with the angles obtainedby HiSCORE (or with generated in ToyMC angles)
void Recon(const char *FileName, int &counter, int draw = 1)
{
    vector<Event> evts;
    ReadTreeOfEvts(FileName, evts);
    int counterOfEvts = 0, counterOfEvts_chi2 = 0;
    for(int ievt = 0; ievt < (int) evts.size(); ievt++)
    {
        if(evts[ievt].NStations() >= N_STN_MIN)
        {
            counterOfEvts++;
            evt.Reset();
            evt = evts[ievt];
            for(int u = 0; u < evt.NStations(); u++) evt.MuStn(u).SetT(evt.MuStn(u).T());
            FillAnglHists(evt, counterOfEvts_chi2);
        }
    }
    cout << "Evts num before chi2 " << counterOfEvts << ", after chi2 " << counterOfEvts_chi2 << endl;
    if(draw) DrawHists();
    counter++;
    cout << "File number = " << counter << endl;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int m()
{
    c52->Divide(2,1);
    int counter = 0;
    Recon("to_analyse/10ns1.root", counter, 0);
    Recon("to_analyse/10ns11.root", counter, 0);
    Recon("to_analyse/10ns12.root", counter, 0);
    Recon("to_analyse/10ns13.root", counter, 0);
    Recon("to_analyse/10ns14.root", counter, 0);
    Recon("to_analyse/10ns15.root", counter, 0);
    Recon("to_analyse/10ns16.root", counter, 0);
    Recon("to_analyse/10ns17.root", counter, 1);
    return 0;
}