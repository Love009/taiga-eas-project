#include "common.h"
using namespace std;
#include <stdlib.h>
#include <bitset>
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "AutoDict_Event_cxx.so"
#else
#include "Tunka.h"
#endif

static const int    N_STN_TOT =  19, N_STN_MIN = 4;
int evt_number;
vector<Event> evts;

vector<float> MuonX;
vector<float> MuonY;
vector<float> MuonZ;
vector<float> MuonTerr;
//============================================================
//Initialize station coordinates and time errors (vectors MuonX, Y, Z, Terr)
// par = 1 for onground stations; par = 2 for underground
void Initialize(int par = 1)
{
    double x, y, z, t, sigma;
    string  title, rest;
    int i = 1;
    ifstream in("coordinatesWithTimeErr.txt");    
    getline(in, title);
    while (1) {
        in >> rest >> x >> y >> z >> sigma;// sigma = t_error
        if (!in.good()) break;
        if(i % 2 == par % 2){//2%2 for getting underground counters coordinstes, 1%2 for onground
            MuonX.push_back(x);
            MuonY.push_back(y);
            MuonZ.push_back(z);
            MuonTerr.push_back(sigma);
        }
        i++;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DrawCol()
{
   Int_t i,n;
   Double_t x, y;
   double t[19], t_[2][19];
   for(i = 0; i < N_STN_TOT; i++) t[i] = -1;
   
  // ReadTreeOfEvts("analysis_results/Gr+Hisc271119_Generation2.root", evts);
   evt_number = 0;
   for(int k = 0; k < 13; k++){
   do{
       //cout << evts[evt_number].PhiGen() << "\t" <<evts[evt_number].ThetaGen()<< endl;
       //if(evts[evt_number].NStations() == N_STN_MIN ) flag++;
       evt_number++;           
    }while(evts[evt_number].NStations() < N_STN_MIN /*|| flag < 10*/);
   }
    cout<<"\t\tevt_number = " << evt_number << "\n\n";
   for(i = 0; i < evts[evt_number].NStations(); i++){
       t[evts[evt_number].MuStn(i).ID() - 1] = evts[evt_number].MuStn(i).T();
   }

   TGraph *g = (TGraph*)gPad->GetListOfPrimitives()->FindObject("Graph");
   n = g->GetN();
   TMarker *m;
   int l = 0;
   for (i = 0; i < n; i++) {       
      g->GetPoint(i,x,y);
      if(t[i] == -1){
          m = new TMarker(x,y,25);
          m-> SetMarkerSize(5);
          m->Paint();
      }
      else{
          t_[0][l] = i;          
          t_[1][l] = t[i];
          l++;
      }
   }   
   for(i = 0; i < l; i++)   cout << t_[1][i] << "\t\t" << t_[0][i] << endl;
    cout<<"\n after:"<<endl;
    double temp1;
    int temp2;//bubble sort
    for (int i = 0; i < l - 1; i++) {
        for (int j = 0; j < l - i - 1; j++) {
            if (t_[1][j] > t_[1][j + 1]) {
                temp1 = t_[1][j];
                temp2 = t_[0][j];
                t_[1][j] = t_[1][j + 1];
                t_[1][j + 1] = temp1;
                t_[0][j] = t_[0][j + 1];
                t_[0][j + 1] = temp2;
            }
        }
    }
    for(i = 0; i < l; i++)   cout << t_[1][i] << "\t\t" << t_[0][i] << endl;
   for(i = 0; i < l; i++)
   {
       g -> GetPoint(t_[0][i], x, y);
       m = new TMarker(x,y,21);
       int color = i*3+57;
       m->SetMarkerColor(color);
       m -> SetMarkerSize(5);
       m->Paint();  
   }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int DrawPicture(){    
    Initialize();   
    ReadTreeOfEvts("ToyMC.root", evts);//analysis_results/Gr+Hisc271119_Generation2.root
    float r[4][19];
    for(int i = 0; i < N_STN_TOT; i++)
    { //for the first generated event
        r[0][i] = MuonX[i];
        r[1][i] = MuonY[i];
        r[2][i] = MuonZ[i];
        r[3][i] = MuonTerr[i];
    }       
    TGraph *gr = new TGraph(19, r[0], r[1]);
    TExec *ex = new TExec("ex","DrawCol();");
    gr->GetListOfFunctions()->Add(ex);
    gr->GetXaxis()->SetLimits(-600,600);
    gr->GetYaxis()->SetRangeUser(-600,600);
    
    gr->GetXaxis()->SetTitle("x, m");
    gStyle->SetTitleY(0.96);
    //gStyle->SetLabelOffset (1);
    gr->GetYaxis()->SetTitleOffset(1.3);
    gr->GetYaxis()->SetTitle("y, m");
    TCanvas *c = new TCanvas("c", "", 1000, 1000);
    gr->SetTitle("Simulated event");//Simulated event
    gr->Draw("AP");
    
    evt_number = 0;
     for(int k = 0; k < 13; k++){
   do{
       //cout << evts[evt_number].PhiGen() << "\t" <<evts[evt_number].ThetaGen()<< endl;
       //if(evts[evt_number].NStations() == N_STN_MIN ) flag++;
       evt_number++;           
    }while(evts[evt_number].NStations() < N_STN_MIN /*|| flag < 10*/);
     }
    cout << "\t\tevt_number = " << evt_number<< "\n\n";
    cout << "phi, theta   " << evts[evt_number].PhiGen() << "\t" <<evts[evt_number].ThetaGen()<< endl;
    double xx=sin(TMath::Pi()/2. + evts[evt_number].ThetaGen())*cos(TMath::Pi()/2. +evts[evt_number].PhiGen());//for real data +pi/2
    double yy=sin(TMath::Pi()/2. + evts[evt_number].ThetaGen())*sin(TMath::Pi()/2. + evts[evt_number].PhiGen());
    double scale = 620.;
    double shift = -50, shiftY=420;
    //double scale = 440.;
    //double shift = 120, shiftY=380;
    TLine *line = new TLine(-scale*xx - shift, -scale*yy + shiftY, scale*xx - shift, scale*yy + shiftY);
    line->SetLineColor(kRed+1);
    line->SetLineWidth(3);
    c->cd();
    c->Range(0.,0.,1.,1.);
    line->Draw();
    //line->DrawLine(0.,0.1,2,2);
    double xx2=sin(evts[evt_number].ThetaGen())*cos(/*TMath::Pi()/2.*3 */+evts[evt_number].PhiGen());
    double yy2=sin(evts[evt_number].ThetaGen())*sin(/*TMath::Pi()/2.*3 */+evts[evt_number].PhiGen());
    scale -=50;
    double scale2 = 150.;
    //double scale2 = 300.;
    TArrow *ar1 = new TArrow( - shift, 0+ shiftY, scale2*xx2 - shift, scale2*yy2+ shiftY, 0.03, "|>");
    //scale2 -= 300;
    TArrow *ar2 = new TArrow(scale*xx - shift, scale*yy+ shiftY, scale*xx+scale2*xx2 - shift, scale*yy+ shiftY+scale2*yy2, 0.03, "|>");//for real data delete +pi/2
    TArrow *ar3 = new TArrow(-scale*xx - shift, -scale*yy+ shiftY, -scale*xx+scale2*xx2 - shift, -scale*yy+ shiftY+scale2*yy2, 0.03, "|>");
    c->cd();
    c->Range(0.,0.,1.,1.);
    ar1->SetLineColor(kRed+1);
    ar1->SetLineWidth(3);
    ar1->SetFillColor(kRed + 1);
    ar1->Draw();
    ar2->SetLineColor(kRed+1);
    ar2->SetLineWidth(3);
    ar2->SetFillColor(kRed + 1);
    ar2->Draw();
    ar3->SetLineColor(kRed+1);
    ar3->SetLineWidth(3);
    ar3->SetFillColor(kRed + 1);
    ar3->Draw();
    
    return 0;
}
