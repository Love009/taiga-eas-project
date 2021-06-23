#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;
#include "TMath.h"

const char *FILE_TO_RECON = "analysis_results/Gr+Hisc261119.root";
static const double SPH_FRONT_RADIUS = 10000.; //meters
bool DEBUG = false;
static const int    N_STN_MIN = 6;
static const int    N_STN_TOT =  19;
static const int    MATCH_HiSCORE = 1;
static const int    SYSTEMATIC = 0;//take in account
double dt_syst[N_STN_TOT];
int matchH = 0;


Event evt;

TH1D *h0   = new TH1D("h0",   "Hist to fit", 20,  0,  20);
TH2D *h1   = new TH2D("h1", "dT_exp-dT_theory [ns], st_i-st_1", 19, 1, 19, 2000, -1000, 1000);
TH2D *phiR   = new TH2D("phiR", "phiR", 1000., 0., 50100., 360., 0., 360.);
TH2D *theR   = new TH2D("theR", "theR", 1000., 0., 50100., 90., 0., 90.);
TH1D *chi2Comparison   = new TH1D("chi2Comparison", "chi2_plane - chi2_sphere", 400., -500., 500.);

TH1D *forTimeResolution   = new TH1D("forTimeResolution", "RMS(t(i)exp-t(i)theory), [ns]", 200., -10., 100.);
TH1D *psi = new TH1D("psi", "psi angle, degrees",   7*10, -2, 15); //angle between reconstructed EAS vector and generated (or reconstructed by HISCORE) one
TH1D *psi_ = new TH1D("psi_", "1 - cos(psi angle)^2",   70*4, 0, 0.067668);
TH1D *dphi = new TH1D("dphi", "delta phi, degrees",   40*2, -20, 20);
TH1D *dthe = new TH1D("dthe", "delta theta, degrees", 40*2, -20, 20);
TH1D *chi2 = new TH1D("chi2", "chi2", 1100, -10, 100);
TH1D *nTrigSts = new TH1D("nTrigSts",   "Number of stations triggered in event", (N_STN_TOT ),  1,  N_STN_TOT);
TH1D *thetaDistr = new TH1D("thetaDistr", "theta, degrees", 190,  -5,  180);
TH1D *phiDistr = new TH1D("phiDistr", "phi, degrees", 370,  -5,  365);
TH1D *thetaDistrRec = new TH1D("thetaDistrRec", "theta, degrees", 190,  -5,  180);
TH1D *phiDistrRec = new TH1D("phiDistrRec", "phi, degrees", 370 / 2.,  -5,  365);
TH1D *cos_the_distrHISC = new TH1D("cos_the_distrHISC", "cos(theta)_HiSCORE", 100., 0., 1.1);
TH1D *cos_the_distrGr = new TH1D("cos_the_distrGr", "cos(theta)_Grande", 100., 0., 1.1);
TH1D *radiusRec = new TH1D("radiusRec","radiusRec", 1000., 0., 50100.);

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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Fit function
Double_t SphericalFrontFit(Double_t *x, Double_t *par)
{
    TRandom2 *randX = new TRandom2();
    TRandom2 *randY = new TRandom2();
    Float_t xx = x[0];
    int ist = int(xx);

    Double_t the = par[0];
    Double_t phi = par[1];
    Double_t t0  = par[2];
    Double_t R0 = par[3];//SPH_FRONT_RADIUS;

    Double_t n1 = sin(the)*cos(phi);
    Double_t n2 = sin(the)*sin(phi);
    Double_t n3 = cos(the);

    double dr[3];
    dr[0] = R0 * n1 + (evt.MuStn(ist).X() - evt.EAScentersGen()[0] - 30. * randX->Gaus(0.,1.));
    dr[1] = R0 * n2 + (evt.MuStn(ist).Y() - evt.EAScentersGen()[1] - 30 * randY->Gaus(0.,1.));
    dr[2] = R0 * n3;
    Double_t val = t0 + sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]) / 0.3 - R0 / 0.3;
    return val;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Fit function 2 (bad polynomial fit)
Double_t FitPsi(Double_t *x, Double_t *par)
{
    Float_t xx = x[0];
    Double_t val = (par[0]+par[1]*xx+par[2]*pow(xx,2)
            +par[3]*pow(xx,3)+par[4]*pow(xx,4)
            +par[5]*pow(xx,5)+par[6]*pow(xx,6)
            +par[7]*pow(xx,7)+par[8]*pow(xx,8)
            +par[9]*pow(xx,9))*xx/par[10];

    return val;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create and fit histogram
void FindDirection(Event& evnt)
{
    double xl = 0.0;
    double xr = evnt.NStations();//number of points of histogram

    TF1  *f0 = new TF1("Plane", PlaneFit, xl, xr, 3);
    h0 -> Reset();
    double T_mean = 0;
    for (int i=0; i < evnt.NStations(); i++) { // Filling the histogram to fit
        double t = evnt.MuStn(i).T();//-evnt.MuStn(0).T();
        T_mean += t;
        double e = evnt.MuStn(i).Terr();
        h0 -> SetBinContent(i+1,t);
        h0 -> SetBinError(i+1, e);
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
            x_max = evnt.MuStn(i).X();
            y_max = evnt.MuStn(i).Y();
        }
    }
    if(DEBUG)cout << "phi_0 = " << 180. / TMath::Pi() * TMath::ATan((y_max-y_min)/(x_max-x_min)) << "(x,y)_min " << x_min << "  " << y_min<< "(x,y)_max " << x_max << " " <<y_max << " T_min, T_max = " << T_min << "   " << T_max <<endl;

    evt = evnt;
    f0->SetParameters(20./180. * TMath::Pi(), TMath::ATan((y_max-y_min)/(x_max-x_min)), T_mean);
    f0->SetParLimits(0, -TMath::Pi() / 2., TMath::Pi() / 2.);
    f0->SetParLimits(1, -TMath::Pi() * 2, TMath::Pi() * 2);
    f0->SetParLimits(2,  -20000, 20000);
   // f0->SetParLimits(2,  0., 1000000000);
    h0->Fit("Plane", "Q0", "", xl, xr); // "Q" - no printing, "P" - Pearson chi-square method
    double the0 = f0 -> GetParameter(0);
    double phi0 = f0 -> GetParameter(1);
    double t00  = f0 -> GetParameter(2);
    double chi0 = f0 -> GetChisquare();

    TF1  *f1 = new TF1("SphericalFront", SphericalFrontFit, xl, xr, 4);
    f1->SetParameters(the0, phi0, t00, 10000.);
    f1->FixParameter(0, the0);
    f1->FixParameter(1, phi0);
    f1->FixParameter(2, t00);
    f1->SetParLimits(0, -TMath::Pi() / 2., TMath::Pi() / 2.);
    f1->SetParLimits(1, -TMath::Pi() * 2, TMath::Pi() * 2);
    f1->SetParLimits(2,  -200000, 200000);
    //f1->SetParLimits(2,  0., 1000000000);
    f1->SetParLimits(3, 0., 30000.);
    h0->Fit("SphericalFront", "Q0", "", xl, xr);
    double the1 = f1 -> GetParameter(0);
    double phi1 = f1 -> GetParameter(1);
    double t01  = f1 -> GetParameter(2);
    double R01  = f1 -> GetParameter(3);

    TF1  *f2 = new TF1("SphericalFront_", SphericalFrontFit, xl, xr, 4);
   // f2->SetParameters(20./180. * TMath::Pi(), TMath::ATan((y_max-y_min)/(x_max-x_min)), T_mean, 3000.);
    f2->SetParameters(the1, phi1, t01, R01);
    f2->SetParLimits(0, -TMath::Pi() / 2., TMath::Pi() / 2.);
    f2->SetParLimits(1, -TMath::Pi() * 2, TMath::Pi() * 2);
    f2->SetParLimits(2,  -200000, 200000);
    //f2->SetParLimits(2,  0., 1000000000);
    f2->SetParLimits(3, 0., 30000.);

//Fitting the histogram
    h0->Fit("SphericalFront_", "Q0", "", xl, xr); // "Q" - no printing, "P" - Pearson chi-square method
    double the = f2 -> GetParameter(0);
    double phi = f2 -> GetParameter(1);

    if(phi < 0. && the > 0.) phi += 2 * TMath::Pi();
    if(the < 0.) {the = -the; phi += TMath::Pi();}
    if(phi < 0.) phi += 2 * TMath::Pi();
    if(phi > TMath::Pi() * 2) phi -= 2 * TMath::Pi();
    double t0  = f2 -> GetParameter(2);
    double chi = f2 -> GetChisquare(); // chi = hi_squared

    chi2Comparison -> Fill(chi0 / (evnt.NStations() - 3.) - chi / (evnt.NStations() - 4.));

    double prp = f2 -> GetProb();
    int    ndf = f2 -> GetNDF();
    double R0 = f2 -> GetParameter(3);
    if(evt.NStations() >= N_STN_MIN && chi < 10.) {radiusRec-> Fill(R0);}
    phiR -> Fill(R0, phi * 57.3);
    theR -> Fill(R0, the * 57.3);
    //else if(evt.NStations() >= N_STN_MIN ) cout << chi << endl;
    //cout << "Chi2/d" << chi << "\t" << evnt.NStations() << endl;
//Writing results
    evnt.SetThetaRec(the);
    evnt.SetPhiRec(phi);
    if(evnt.NStations() > 4) evnt.SetChi2(chi / (evnt.NStations() - 4.));

    int *TimeCust;
    TimeCust = evnt.MuStn(stIDmin).Trow();
    TimeCust[3] = (int) t0 - evnt.MuStn(stIDmin).Trow()[3];
   // cout <<"TimeCust[3] = " << TimeCust[3] << endl;
    evnt.SetEvtTime(TimeCust);

    if(evnt.NStations() == N_STN_MIN)
    {
        double tk = 0.;
        for(int k = 0; k < evnt.NStations(); k++)
        {
            double x, y, z;
            x = evnt.MuStn(k).X();
            y = evnt.MuStn(k).Y();
            z = evnt.MuStn(k).Z();
            double n1 = cos(evnt.PhiRec())*sin(evnt.ThetaRec());
            double n2 = sin(evnt.PhiRec())*sin(evnt.ThetaRec());
            double n3 = cos(evnt.ThetaRec());
            tk += pow(evnt.MuStn(k).T() - (n1*x+n2*y+n3*z)/0.3 - t0,2);
        }
        forTimeResolution -> Fill(sqrt(tk));
    }
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
    auto c000 = new TCanvas("c000", "");
    chi2Comparison -> Draw();
    auto c00 = new TCanvas("c00","");
    radiusRec -> Draw();
    auto c188 = new TCanvas("c188", "");
    psi->Draw();
    /*auto c1 = new TCanvas("c1", "");
    c1->Divide(1,2);
    c1->cd(1);
    //gr -> Draw("AP*");
    psi_->Draw();
    gStyle->SetOptStat(111);//output without RMS
    c1->cd(2);
    //auto c1 = new TCanvas("c1", "");
    //c1->SetLogy();
    int interval = 8;
   // TF1 *fit_psi = new TF1("fit_psi", "pol9", 0, interval);
  //  psi->Fit("fit_psi","","", 0, interval);
    psi->Draw();
    auto c188 = new TCanvas("c188", "");
    psi->Draw();*/

    auto c2 = new TCanvas("c2", "");
    c2->Divide(2, 1);
    c2->cd(1);
    //gPad->SetLogy();
    dphi->Draw();
    c2->cd(2);
    dthe->Draw();

    auto c218 = new TCanvas("c218","");
    TF1 *exp_cust = new TF1("exp_cust", "[1] * exp(-x/[0])", 0., 20.);
    exp_cust -> SetParLimits(0, 0.001, 30.);
    double norm10 = exp_cust -> GetMaximum();
    exp_cust -> SetParLimits(1, 0., 1000000);
    //exp_cust ->FixParameter(1, norm10);
   // exp_cust ->FixParameter(0, 2.);
    chi2 -> Fit("exp_cust", "","",0., 20.);
    chi2 -> Draw("E0");
// For theta real distr drawing
    auto c11 = new TCanvas("c11", "");
    cos_the_distrHISC -> Scale(1. / cos_the_distrHISC -> GetEntries());
    cos_the_distrHISC -> GetXaxis() -> SetTitle("cos(theta)");
    double norm11 = cos_the_distrHISC -> GetMaximum();
    TF1 *pol_cos = new TF1("pol_cos", "[1] * pow(x,[0])", 0.57, 1.);
    pol_cos -> SetParLimits(0, 2.,20.);
    pol_cos ->FixParameter(1, norm11);
    //cos_the_distrHISC -> Fit("pol_cos", "","",0.57,1.);
   // gStyle->SetOptFit(1111);
    cos_the_distrHISC -> Draw();

    auto c12 = new TCanvas("c12", "");
    cos_the_distrGr -> Scale(1. / cos_the_distrGr -> GetEntries());
    cos_the_distrGr -> GetXaxis() -> SetTitle("cos(theta)");
    double norm12 = cos_the_distrGr -> GetMaximum();
    TF1 *pol_cos2 = new TF1("pol_cos2", "[1] * pow(x,[0]) + [2]", 0., 1.);
    pol_cos2 -> SetParLimits(0, -2.,20.);
    pol_cos2 ->FixParameter(1, norm12);
    pol_cos2 -> SetParLimits(2, 0.,0.5);
    cos_the_distrGr -> Fit("pol_cos2", "","",0.,1.);
    gStyle->SetOptFit(1111);
    cos_the_distrGr -> Draw();
//For real and gen phi and theta distr drawing
    auto c15 = new TCanvas("c15", "");
    //phiDistr -> Draw();
    //phiDistrRec -> SetLineColor(kBlack);
    TF1 *constant = new TF1("constant", "[0]", 0., 360.);
    constant-> SetParLimits(0, 0.,200000.);
   // phiDistrRec -> Fit("constant", "","",0.,360.);
    phiDistrRec -> Draw("E0");
    auto c16 = new TCanvas("c16", "");
    thetaDistr -> Draw();
    thetaDistrRec -> SetLineColor(kRed);
    thetaDistrRec -> Draw("same");
    auto c5 = new TCanvas("c5", "");
    h1->SetMarkerStyle(kFullTriangleDown);
    h1->Draw();
    auto c100 = new TCanvas("c100","");
    forTimeResolution -> Draw();

    auto c115 = new TCanvas("c115","");
    phiR->Draw();
    auto c120 = new TCanvas("c120","");
    theR->Draw();
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CheckSpecificInterval(Event &evt)
{
    double center = 45. / 180. * TMath::Pi();
    double interval = 10. / 180. * TMath::Pi() / 2.;
    if(evt.ThetaRec() <= center + interval && evt.ThetaRec() >= center - interval) return 1;
    else return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FillAnglHists(Event &evt, int &counterOfEvts_chi2)
{
    int toOptimize = 1;
    do{
        FindDirection(evt);
       // cout << 180. / TMath::Pi()*evt.PhiGen() - 180. / TMath::Pi()*evt.PhiRec() <<"\t" << 180. / TMath::Pi()*evt.ThetaGen() - 180. / TMath::Pi()*evt.ThetaRec() << endl;
        double d_psi = sin(evt.ThetaRec())*cos(evt.PhiRec())*sin(evt.ThetaGen())*cos(evt.PhiGen())
                       + sin(evt.ThetaRec())*sin(evt.PhiRec())*sin(evt.ThetaGen())*sin(evt.PhiGen())
                       + cos(evt.ThetaRec())*cos(evt.ThetaGen());

        if(evt.Chi2() < 10)
        {
            toOptimize = 0;
            counterOfEvts_chi2++;
            if (DEBUG)    cout  << "Reconstructed angle Theta = " << evt.ThetaRec() * 180. / TMath::Pi() <<  " : " << evt.ThetaGen() * 180. / TMath::Pi() << "\n"  << "Reconstructed angle Phi   = "   << evt.PhiRec() * 180. / TMath::Pi()  <<  " : " <<evt.PhiGen() * 180. / TMath::Pi()  << "\n"  << "Chi2 of the fit           = "            << evt.Chi2()     <<"\n" << endl;

            phiDistrRec->Fill(180. / TMath::Pi() * evt.PhiRec());

            if(evt.MatchHiSC() == MATCH_HiSCORE) {
                matchH++;
                double temp0 = TMath::ACos(d_psi);
                if(temp0 > TMath::Pi() / 2.) temp0 = TMath::Pi() - temp0;

                /*if(CheckSpecificInterval(evt))*/psi -> Fill(180. / TMath::Pi() * temp0);
                psi_ -> Fill(1 - pow(TMath::Cos(temp0),2));
                thetaDistr->Fill(180. / TMath::Pi()*evt.ThetaGen());
                thetaDistrRec->Fill(180. / TMath::Pi()*evt.ThetaRec());
                dthe -> Fill(180. / TMath::Pi() * (evt.ThetaRec() - evt.ThetaGen()));
                dphi -> Fill(180. / TMath::Pi() * (evt.PhiRec() - evt.PhiGen())*TMath::Sin(evt.ThetaRec()));
                phiDistr->Fill(180. / TMath::Pi() * evt.PhiGen());
               // phiDistrRec->Fill(180. / TMath::Pi() * evt.PhiRec());
                cos_the_distrHISC -> Fill(TMath::Cos(evt.ThetaGen()));
                cos_the_distrGr -> Fill(TMath::Cos(evt.ThetaRec()));
            }

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
    if(evt.NStations() >= 4) chi2 -> Fill(evt.Chi2());
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FillTimeHists(Event &evt)
{
    for(int j = 0; j < evt.NStations(); j++)
    {
        for(int k = 0; k < evt.NStations(); k++)
        {
            if(k != j && evt.NStations() >= 3)
            {
                double dT = FindDT(evt, k, j);
                if(evt.MuStn(j).ID() == 1) h1 -> Fill(evt.MuStn(k).ID(), dT);
            }
        }
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static const double c = 0.3; //speed of light, ns
static const int    N_STN_MIN_SYS = 2;
vector<double> MuonX;
vector<double> MuonY;
vector<double> MuonZ;
const char*    FILE_NAME    = "ToyMC.root";
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
    //for(int u = 0; u < N_STN_TOT; u++) M_tilda_tot(0, u) = 0.;
   // M_tilda_tot(0, 18) = 1.;
    res(0, 0) = 0.;

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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CheckHiSCAperture(Event &evt)
{
    double HiSC_the = 20. / (180. / TMath::Pi());
    double HiSC_view_size = 30. / (180. / TMath::Pi());
    double the = evt.ThetaGen() , phi = evt.PhiGen();
    if(TMath::ACos(sin(HiSC_the) * sin(the) * cos(phi) + cos(HiSC_the) * cos(the)) <= HiSC_view_size) return 1;
    else return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CheckCenter600(Event &evt)
{
    if(sqrt(evt.EAScentersGen()[0] * evt.EAScentersGen()[0] + evt.EAScentersGen()[1] * evt.EAScentersGen()[1]) <= 600.) return 1;
    else return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Reconstruct phi and theta angles of EASes and compare with the angles obtainedby HiSCORE (or with generated in ToyMC angles)
void Recon(const char *FileName = FILE_TO_RECON, int draw = 1)
{
    double *temp;
    if(SYSTEMATIC) {
        temp = FindTSystematic(FileName);
        for (int u = 0; u < N_STN_TOT; u++) dt_syst[u] = *(temp + u);
    }
    else
    {
        for (int u = 0; u < N_STN_TOT; u++) dt_syst[u] = 0.;
    }
    //for(int u = 0; u < N_STN_TOT; u++) dt_syst[u] = 0.;

    vector<Event> evts;
    ReadTreeOfEvts(FileName, evts);
    int counterOfEvts = 0, counterOfEvts_chi2 = 0;
    for(int ievt = 0; ievt < (int) evts.size(); ievt++)
    {
        evts[ievt].SetPhiRec(-100.);
        if(evts[ievt].NStations() > 2) nTrigSts->Fill(evts[ievt].NStations());
        int flag600 = 1;//CheckCenter600(evts[ievt]);
        int flagHiSCAperture = 1;//CheckHiSCAperture(evts[ievt]);
        if(evts[ievt].NStations() >= N_STN_MIN && flag600 == 1 && flagHiSCAperture == 1)
        {
            counterOfEvts++;
            for(int u = 0; u < evts[ievt].NStations(); u++) evts[ievt].MuStn(u).SetT(  evts[ievt].MuStn(u).T() - dt_syst[evts[ievt].MuStn(u).ID() - 1]);
            FillAnglHists(evts[ievt], counterOfEvts_chi2);
            FillTimeHists(evts[ievt]);
        }
    }
    cout << "Evts num before chi2 " << counterOfEvts << ", after chi2 " << counterOfEvts_chi2 << endl;
    cout <<  matchH << " MATCHED evts" << endl;
    WriteToFile(evts, FileName);
    if(draw) DrawHists();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Some custom functions for quick analysis:
int main()
{
    //Recon( "ToyMC.root");
    Recon("analysis_results/Gr+Hisc261119_Generation3.root");
    Recon("analysis_results/Gr+Hisc271119_Generation3.root");
    TH1D *gen_phi = (TH1D*) thetaDistr->Clone();
    TH1D *gen_the = (TH1D*) phiDistr->Clone();
    TH1D *rec_phi = (TH1D*) thetaDistrRec->Clone();
    TH1D *rec_the = (TH1D*) phiDistrRec->Clone();
    Double_t norm = gen_phi->GetEntries();
    auto c7 = new TCanvas("c7","");
    rec_phi->Draw("");
    gen_phi->SetLineColor(kRed);
    gen_phi->Draw("HIST SAME");
    auto c8 = new TCanvas("c8","");
    rec_the->Draw("");
    gen_the->SetLineColor(kRed);
    gen_the->Draw("HIST SAME");
    return 0;
}

int main3()
{
    Recon("analysis_results/for_EAS_center/Gr+Hisc241119.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc241219.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc291219.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc020120.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc221119.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc261119.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc271119.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc301219.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc150220.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc160220.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc200220.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc250220.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc270220.root", 0);
    Recon("analysis_results/for_EAS_center/Gr+Hisc210220.root", 1);
    return 0;
}

