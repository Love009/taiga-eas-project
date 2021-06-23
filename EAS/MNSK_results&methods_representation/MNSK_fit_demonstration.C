#include "common.h"
using namespace std;
#include <stdlib.h>
#include <bitset>
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "AutoDict_Event_cxx.so"
#else
#include "Tunka.h"
#endif

bool DEBUG = false;

Event evt;
static const int    N_STN_TOT =  19;
static const int    N_STN_MIN =   3;
static const int    N_STN_MAX =  19;
static const double P_STN_TRG = 0.3;  

vector<float> MuonX;
vector<float> MuonY;
vector<float> MuonZ;
vector<float> MuonTerr;

//TH1D *h0   = new TH1D("h0",   "Hist to fit", 20,  0,  20);
auto c0 = new TCanvas("c0","");


//TH1D *h0   = new TH1D("h0",   "Hist to fit", 20,  0,  20);
TH2D *h1   = new TH2D("h1", "dT_exp-dT_theory [ns], st_i-st_1", 19, 1, 19, 2000, -1000, 1000);
TH1D *h2   = new TH1D("h2", "dT_exp-dT_theory [ns], 2 neighbour stations", 2000, -1000, 1000);
TH1D *dt_hist   = new TH1D("dt_hist",   "dT_exp - dT_theory [ns]", 1000, -500, 500);
TH1D *dpsi = new TH1D("dpsi", "psi angle, degrees",   120, -10, 50); //angle between reconstructed EAS vector and generated (or reconstructed by HISCORE) one
TH1D *dphi = new TH1D("dphi", "delta phi, degrees",   200, -50, 50);
TH1D *dthe = new TH1D("dthe", "delta theta, degrees", 200, -50, 50);
TH1D *chi2 = new TH1D("chi2", "chi2", 110, -10, 100);

//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
void Initialize()
{
    double x, y, z, t, sigma;
    string  title, rest;
    int i = 0;
    ifstream in("coordinatesWithTimeErr.txt");    
    getline(in, title);
    while (1) {
        in >> rest >> x >> y >> z >> sigma;// sigma = t_error
        if (!in.good()) break;
        if(i % 2 == 0){
            MuonX.push_back(x);
            MuonY.push_back(y);
            MuonZ.push_back(z);
            MuonTerr.push_back(sigma);
        }
        i++;
    }
}

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


void FindDirection(Event& evt)
{
  double xl = 0.0;
  double xr = evt.NStations();//number of points of histogram

  TF1  *f1 = new TF1("Plane", PlaneFit, xl, xr, 2);
  TH1D *h0 = new TH1D("h0", "Hist to fit",   20, 0, 20.);

  //h0 -> Reset();
  for (int i=0; i < evt.NStations(); i++) { // Filling the histogram to fit

    double t = evt.MuStn(i).T();
    double e = evt.MuStn(i).Terr();

    h0 -> SetBinContent(i+1,t);
    h0 -> SetBinError(i+1,e);
    RooFit::SumW2Error(1);
  }

 f1->SetParameters(0.2, 3.0, 1000.);
  f1->SetParLimits(0, -1.58, 1.58);
  f1->SetParLimits(1, -6.3, 6.3); //TMath::Pi() * 2);
  f1->SetParLimits(2,  -20000, 20000);
  f1->SetLineStyle(3);

  //h0->Fit("Plane", "NQP", "", xl, xr);          // Fitting the histogram 
  h0->Fit("Plane", "Q", "", xl, xr); // "Q" - no printing, "P" - Pearson chi-square method
  double the = f1 -> GetParameter(0);
  double phi = f1 -> GetParameter(1);
  double chi = f1 -> GetChisquare(); // chi = hi_squared
  double prp = f1 -> GetProb();
  int    ndf = f1 -> GetNDF();
  
  evt.SetThetaRec(the);
  evt.SetPhiRec(phi);
  evt.SetChi2(chi);
  
  if (DEBUG) cout << " Number of stations: " << evt.NStations() << endl;
  if (DEBUG) for (int i=0; i<evt.NStations(); i++) {
      cout << "           station: " << evt.MuStn(i).ID() << "  " << 
evt.MuStn(i).T() << endl;
    }
  f1->SetLineWidth(1);
  h0->Draw();
  f1 -> Draw("SAMEalf");  
  c0->Modified();
  RooFit::SumW2Error(1);
  c0->Update();  
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Read vector of Events from a root-file and put it in evts
void ReadTreeOfEvts(const char *FileName, vector<Event> &evts)
{
    TFile *f = new TFile(FileName);
    TTree *t4 = (TTree*)f->Get("Muon");
    Event *event = new Event();
    TBranch *br = t4->GetBranch("events");
    br ->SetAddress(&event);
    long int nEvTot = t4 -> GetEntries();
    for(long int i = 0; i < nEvTot; i++)
    {
       br ->GetEntry(i);
       evts.push_back(*event);
       event->Reset();
    }
    cout << nEvTot << " events" << endl;
    f->Close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Reconstruct phi and theta angles of EASes and compare with the angles obtainedby HiSCORE (or with generated in ToyMC angles)
void Recon(const char *FileName = "ToyMC.root")
{
    vector<Event> evts, evtsNew;
    ReadTreeOfEvts(FileName, evts);    

    for(int ievt = 0; ievt < (int) evts.size(); ievt++)
    {
        if(evts[ievt].NStations() >= N_STN_MIN)
        {
            evt.Reset();
            evt = evts[ievt];
            FindDirection(evt);          
            
                
            double d_psi = sin(evt.ThetaRec())*cos(evt.PhiRec())*sin(evt.ThetaGen())*cos(evt.PhiGen()) 
                                    + sin(evt.ThetaRec())*sin(evt.PhiRec())*sin(evt.ThetaGen())*sin(evt.PhiGen()) 
                                    + cos(evt.ThetaRec())*cos(evt.ThetaGen()); 
            double temp = evt.PhiRec() - evt.PhiGen();
            if(temp >= TMath::Pi()) temp = 2 *TMath::Pi() - temp;
            if(temp <= -TMath::Pi()) temp = 2 * TMath::Pi() + temp;
            if(evt.Chi2() < 10)
            {
                if (DEBUG)    cout  << "Reconstructed angle Theta = " << evt.ThetaRec() * 57.3 <<  " : " << evt.ThetaGen() * 57.3 << "\n"  << "Reconstructed angle Phi   = "   << evt.PhiRec() * 57.3  <<  " : " <<evt.PhiGen() * 57.3  << "\n"  << "Chi2 of the fit           = "            << evt.Chi2()     <<"\n" << endl;
                dpsi -> Fill(57.3 * TMath::ACos(d_psi)); 
                double ttemp = evt.ThetaRec();
                if(evt.ThetaRec() < 0) ttemp = ttemp * (-1);
                dthe -> Fill(57.3 * (ttemp - evt.ThetaGen()));
                dphi -> Fill(57.3 * temp*TMath::Sin(evt.ThetaRec()));
                
            }
            if(evt.NStations() > 3) chi2 -> Fill(evt.Chi2());
            
//Histograms filling for representation of the reconstruction results in terms dT
            for(int j = 1; j < evt.NStations(); j++)
            {
                int min_j = 0;
                double dt_exp = double(evt.MuStn(j).T()- evt.MuStn(min_j).T());
                double x, y, z;
                x = evt.MuStn(j).X() - evt.MuStn(min_j).X();
                y = evt.MuStn(j).Y() - evt.MuStn(min_j).Y();
                z = evt.MuStn(j).Z() - evt.MuStn(min_j).Z();
                double n1 = cos(evt.PhiRec())*sin(evt.ThetaRec());
                double n2 = sin(evt.PhiRec())*sin(evt.ThetaRec());
                double n3 = cos(evt.ThetaRec());
                double dt_theory = (n1*x+n2*y+n3*z)/0.3;
                double dT = dt_exp - dt_theory;
                if(evt.NStations() >= 3)    dt_hist->Fill(dT);           
            }  
            for(int j = 0; j < evt.NStations(); j++)
            {            
                for(int k = 0; k <  evt.NStations(); k++)
                {
                    if(k != j)
                    {
                        double dt_exp = double(evt.MuStn(k).T()- evt.MuStn(j).T());
                        double x, y, z;
                        x = evt.MuStn(k).X() - evt.MuStn(j).X();
                        y = evt.MuStn(k).Y() - evt.MuStn(j).Y();
                        z = evt.MuStn(k).Z() - evt.MuStn(j).Z();
                        double n1 = cos(evt.PhiRec())*sin(evt.ThetaRec());
                        double n2 = sin(evt.PhiRec())*sin(evt.ThetaRec());
                        double n3 = cos(evt.ThetaRec());
                        double dt_theory = (n1*x+n2*y+n3*z)/0.3;
                        double dT = dt_exp - dt_theory;
                        if(evt.NStations() >= 3)
                        {
                            dt_hist->Fill(dT);   
                            if(evt.MuStn(j).ID() == 1) h1 -> Fill(evt.MuStn(k).ID(), dT);
                            if(-evt.MuStn(j).ID() + evt.MuStn(k).ID() == 1) h2 -> Fill(dT);
                        }
                    }
                }            
            }
        }
    }
    /*auto c1 = new TCanvas("c1","");
    //c1->SetLogy();
    dpsi->Draw();
    auto c2 = new TCanvas("c2","");
    c2 -> Divide(2,1);
    c2 -> cd(1);
    //gPad->SetLogy();
    dphi->Draw();
    c2 -> cd(2);
    //gPad->SetLogy();
    dthe->Draw();
    auto c3 = new TCanvas("c3","");
    gPad->SetLogy();
    chi2 ->Draw();
    auto c4 = new TCanvas("c4","");
    gPad->SetLogy();
    dt_hist -> Draw();
    auto c5 = new TCanvas("c5","");
    h1->SetMarkerStyle(kFullTriangleDown);
    h1->Draw();  
    auto c6 = new TCanvas("c6","");
    gPad->SetLogy();
    h2->Draw();*/
}

/*void Recon(const char *FileName = "ToyMC.root", int evtType = 1)
{
  
#include "RootStruct.h"

  TFile *f = new TFile(FileName);
  TTree *data = (TTree*)f->Get("Muon");

  data->SetBranchAddress("nstn",    &nstn);
  data->SetBranchAddress("stID",     stID);
  data->SetBranchAddress("stX",       stX);
  data->SetBranchAddress("stY",       stY);
  data->SetBranchAddress("stZ",       stZ);
  data->SetBranchAddress("stT",       stT);
  data->SetBranchAddress("stTerr", stTerr);

  if (evtType==1){
    data->SetBranchAddress("the_gen",    &the_gen);
    data->SetBranchAddress("phi_gen",    &phi_gen);
  }
  
  int Nevent = data->GetEntries();
  cout << "Number of events = " << Nevent << endl;

  for (int ievt=0; ievt<Nevent; ievt++)
    {
      data->GetEntry(ievt);

      evt.Reset();
      evt.SetThetaGen(the_gen);
      evt.SetPhiGen(phi_gen);
      
      for (int istn=0; istn<nstn; istn++) {
	Station stn;
	
	stn.SetID(stID[istn]);
	stn.SetX(stX[istn]);
	stn.SetY(stY[istn]);
	stn.SetZ(stZ[istn]);
	stn.SetT(stT[istn]);
	stn.SetTerr(stTerr[istn]);

	evt.AddStation(stn);
      }

      if (evt.NStations()>N_STN_MIN && evt.NStations()<N_STN_MAX)
	{  int numb = 0;    
	  FindDirection(evt);

	  if (DEBUG) cout << "Reconstructed angle Theta = " << evt.ThetaRec() << " : " << evt.ThetaGen() << endl;
	  if (DEBUG) cout << "Reconstructed angle Phi   = " << evt.PhiRec()   << " : " << evt.PhiGen()   << endl;
	  if (DEBUG) cout << "Chi2 of the fit           = " << evt.Chi2()     << endl;

	  double d_psi = sin(evt.ThetaRec())*cos(evt.PhiRec())*sin(evt.ThetaGen())*cos(evt.PhiGen()) +
	                 sin(evt.ThetaRec())*sin(evt.PhiRec())*sin(evt.ThetaGen())*sin(evt.PhiGen()) +
	                 cos(evt.ThetaRec())*cos(evt.ThetaGen()); 
	  dpsi -> Fill(57.3*TMath::ACos(d_psi));
	  dthe -> Fill(57.3*(evt.ThetaRec() - evt.ThetaGen()));
	  dphi -> Fill(57.3*(evt.PhiRec() - evt.PhiGen()));
	  chi2 -> Fill(evt.Chi2());
	}
    }
  
  f->Close();
  auto c1 = new TCanvas("c1",""); 
  dpsi->Draw();
  auto c2 = new TCanvas("c2",""); 
  c2 -> Divide(2,1);
  c2 -> cd(1);
  dphi->Draw();
  c2 -> cd(2);
  dthe->Draw();
  auto c3 = new TCanvas("c3","");
  chi2 ->Draw();

  
  f->Close();
}*/

