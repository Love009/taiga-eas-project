#include "common.h"
using namespace std;
#include <stdlib.h>
#include <bitset>
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "AutoDict_Event_cxx.so"
#else
#include "TunkaIrkutsk.h"
#endif

static const double P_STN_TRG = 1;
static const double R0 = 120.; //meters
const char*         FILE_NAME    = "ToyMC.root";

vector<double> MuonX[2];
vector<double> MuonY[2];
vector<double> MuonZ[2];
vector<double> MuonTerr[2];
//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
//Initialize station coordinates and time errors (vectors MuonX, Y, Z, Terr)
// par = 1 for onground stations; par = 2 for underground
void Initialize()
{
    double x, y, z, t, sigma;
    string  title, rest;
    int i = 1;
    ifstream in("coordinatesWithTimeErr.txt");    
    getline(in, title);
    while (1) {
        in >> rest >> x >> y >> z >> sigma;// sigma = t_error
        if (!in.good()) break;
        if(i % 2 == 1){//2%2 for getting underground counters coordinstes, 1%2 for onground
            MuonX[0].push_back(x);
            MuonY[0].push_back(y);
            MuonZ[0].push_back(z);
            MuonTerr[0].push_back(sigma);
        }
        else{
            MuonX[1].push_back(x);
            MuonY[1].push_back(y);
            MuonZ[1].push_back(z);
            MuonTerr[1].push_back(sigma);
        }
        i++;
    }
}
//===============================================================================================//
//Events generation and writing to the root-file
int Generate(int Ngen = 10000, double syst = 0., const char* fileToSave = FILE_NAME)
{
    double t_syst[19];
    for(int i = 0; i < 19; i++) t_syst[i] = 0.;
    t_syst[0] = syst;
    TH1D *nTrigSts = new TH1D("nTrigSts",   "N triggered stations", (N_STN_TOT + 1),  0,  N_STN_TOT); 
    TH1D *xx = new TH1D("xx",   "X coord", 1200,  -600,  600); 
    TH1D *phi = new TH1D("phi",   "phi", 370,  -5,  365); 
    TH1D *theta = new TH1D("theta",   "theta", 190,  -5,  180); 
    TH1D *prob = new TH1D("prob",   "Prob-ty", 1000,  0,  1); 
    int nstn;
    double the_gen, phi_gen, cos_gen[3], x_coord, y_coord; 
    TFile f(fileToSave,"recreate");
    TTree Muon("Muon","");
    Event *event = new Event();
    Muon.Branch("events", &event);     
    TRandom2 *r2=new TRandom2();
    
    Initialize(); 
    
    for (int i = 0; i < Ngen; i++)
    {
        event->Reset();
        the_gen    =TMath::ASin(gRandom -> Rndm());
        phi_gen    = 2*TMath::Pi()*(gRandom -> Rndm());
        phi->Fill(57.3*phi_gen);
        theta->Fill(57.3*the_gen);
        event->SetPhiHiSC(phi_gen);
        event->SetThetaHiSC(the_gen);
        
        cos_gen[0] = sin(the_gen)*cos(phi_gen);
        cos_gen[1] = sin(the_gen)*sin(phi_gen);
        cos_gen[2] = cos(the_gen);
        x_coord = 1200. * (gRandom -> Rndm()) - 600.;
        xx->Fill(x_coord);
        y_coord = 1200. * (gRandom -> Rndm()) - 600.;

        //Initialize();
        nstn = 0;
        for (int j = 0; j < N_STN_TOT; j++) {
            double R = sqrt((MuonX[0][j] - x_coord)*(MuonX[0][j] - x_coord) + (MuonY[0][j] - y_coord)* (MuonY[0][j] - y_coord));
            prob->Fill(P_STN_TRG * exp(-R/R0));
            if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) {
                Station st;
                double x[2], y[2], z[2];
                x[0] = MuonX[0][j];
                x[1] = MuonX[1][j];
                y[0] = MuonY[0][j];
                y[1] = MuonY[1][j];
                z[0] = MuonZ[0][j];
                z[1] = MuonZ[1][j];
                int Traw[4];
                Traw[3] = (cos_gen[0]*MuonX[0][j]  + cos_gen[1]*MuonY[0][j] + cos_gen[2]*MuonZ[0][j]) / 0.3 + MuonTerr[0][j]*(r2->Gaus(0.,1.)) + t_syst[j];
                
                st.SetID(j + 1);
                st.SetY(y);
                st.SetX(x);
                st.SetZ(z);
                st.SetTraw(Traw);
                event->AddStation(st);
                nstn++;
            }
        }
        if(nstn>2) nTrigSts -> Fill(nstn);
        Muon.Fill();
    }
    f.Write();
    f.Close();
    return 0;
}
