#include "common.h"
using namespace std;
#include <stdlib.h>
#include <bitset>
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "AutoDict_Event_cxx.so"
#else
#include "Tunka.h"
#endif

static const int    N_STN_TOT =  19;
static const double P_STN_TRG = 1;
static const double R0 = 120.; //meters
const char*         FILE_NAME    = "analysis_results/ToyMC.root";
double              sep = 200000.;

vector<float> MuonX;
vector<float> MuonY;
vector<float> MuonZ;
vector<float> MuonTerr;
//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
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
//===============================================================================================//
//Events generation and writing to the root-file
int Generate(int Ngen = 1000, const char* fileToSave = FILE_NAME)
{
    TH1D *phi = new TH1D("phi",   "phi", 370,  -5,  365);
    TH1D *theta = new TH1D("theta",   "theta", 360,  -180,  180);
    /*TH1D *nTrigSts = new TH1D("nTrigSts",   "N triggered stations", (N_STN_TOT + 1),  0,  N_STN_TOT);
    TH1D *xx = new TH1D("xx",   "X coord", 1200,  -600,  600); 

    TH1D *prob = new TH1D("prob",   "Prob-ty", 1000,  0,  1); */
    int nstn;
    double the_gen, phi_gen, cos_gen[3], x_coord, y_coord;
   // TFile f(fileToSave,"recreate");
    TTree Muon("Muon","");
    Event *event = new Event();
    Muon.Branch("events", &event);     
    TRandom2 *r2=new TRandom2();

    Initialize(1);
    for (int i = 0; i < Ngen; i++)
    {
    //    int s = 1;
        event->Reset();
        double moon_phi = 100. / 57.3;
        double moon_the = 35. / 57.3;
        double range_moon = 180. / 57.3;
        double moon_size =0.;// 0.25 / 57.3;
        int flag = 0;
        //do {
            flag = 0;
            the_gen = TMath::ASin(gRandom->Rndm());
            //if(the_gen < 0) the_gen = -the_gen;
           /* double y;
            do
            {
                the_gen = TMath::Pi() / 2. * (gRandom->Rndm());
                //the_gen = (2 * range_moon / 180.) * TMath::Pi() * (gRandom->Rndm()) + (moon_the - range_moon) / 180. * TMath::Pi();
                if(the_gen >= TMath::Pi() / 2.) the_gen -= TMath::Pi() / 2.;
                y = (gRandom->Rndm());
            }while(y > pow(TMath::Cos(the_gen), 8));*/
            phi_gen = (2 * range_moon) * (gRandom->Rndm());// + (moon_phi - range_moon);
           // phi_gen = 2 * TMath::Pi() * (gRandom->Rndm());
           // s = TMath::Sin(the_gen) / abs(TMath::Sin(the_gen));
            if((the_gen - moon_the) * (the_gen - moon_the) + (phi_gen - moon_phi) * (phi_gen - moon_phi) * TMath::Sin(moon_the) * TMath::Sin(moon_the) <= moon_size * moon_size ) flag = 1;
            else flag = 0;
       // }while(flag);

       phi->Fill(57.3*phi_gen);
        theta->Fill(57.3*the_gen);
        event->SetPhiGen(phi_gen);
        event->SetThetaGen(the_gen);
        
        cos_gen[0] = sin(the_gen)*cos(phi_gen);
        cos_gen[1] = sin(the_gen)*sin(phi_gen);
        cos_gen[2] = cos(the_gen);
        x_coord = 1200. * (gRandom -> Rndm()) - 600.;
       // xx->Fill(x_coord);
        y_coord = 1200. * (gRandom -> Rndm()) - 600.;

        //Initialize();
        nstn = 0;
        for (int j = 0; j < N_STN_TOT; j++) {
            double R = sqrt((MuonX[j] - x_coord)*(MuonX[j] - x_coord) + (MuonY[j] - y_coord)* (MuonY[j] - y_coord));
           // prob->Fill(P_STN_TRG * exp(-R/R0));
            if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) {
                Station st;
                st.SetID(j + 1);
                st.SetX(MuonX[j]);
                st.SetY(MuonY[j]);
                st.SetZ(MuonZ[j]);
                st.SetT((cos_gen[0]*MuonX[j]  + cos_gen[1]*MuonY[j] + cos_gen[2]*MuonZ[j]) / 0.3 + /*MuonTerr[j]*/15. *(r2->Gaus(0.,1.)));
                st.SetTerr(MuonTerr[j]);
                event->AddStation(st);
                event->SetMatchHiSC(1);
                nstn++;
            }
        }

        if(nstn > 2)
        {
          //  nTrigSts -> Fill(nstn);phi->Fill(57.3*phi_gen);
           // theta->Fill(57.3*the_gen);
            Muon.Fill();

            if((i + 1) / sep - (int)((i + 1) / sep) == 0 && i != 0)
            {
                cout << i << endl;
                int num = (int)(i / sep);
                string name(fileToSave);
                name = name + to_string(num) + ".root";
                cout << name << endl;
                TFile f(name.c_str(),"recreate");
                Muon.Write();
                f.Close();
                Muon.Reset();
            }
            else if(i == Ngen - 1)
            {
                TFile f(fileToSave,"recreate");
                Muon.Write();
                f.Close();
            }
        }
        else i--;//!!!!!!!!!!!!!!!!!!!!!!!!!only here it is necessary
        //Muon.Fill();
    }

    //f.Write();
    //f.Close();
    //Double_t norm = nTrigSts->GetEntries();
    //nTrigSts->Scale(1/norm);
    
    //nTrigSts -> Draw();
   // auto c2 = new TCanvas("c2","");
    //nTrigSts -> Draw();
   // auto c3 = new TCanvas("c3","");
  //  prob -> Draw();
    auto c4 = new TCanvas("c4","");
    phi -> Draw();
    auto c5 = new TCanvas("c5","");
    theta -> Draw();
    return 0;
}
