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
static const double R0 = 150.; //meters
const char*         FILE_NAME    = "analysis_results/ToyMC.root";
double              sep = 200000.;
double              TIME_ERROR = 5.; //ns

vector<float> MuonX;
vector<float> MuonY;
vector<float> MuonZ;
vector<float> MuonTerr;
int n_bins = 25 * 10;
TH2D *theta_phi_gen = new TH2D("theta_phi_gen", "theta phi gen, degrees", n_bins,  -2.,  2., n_bins,  -2., 2.);
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
    int nstn;
    double the_gen, phi_gen, cos_gen[3], x_coord, y_coord;
    TTree Muon("Muon","");
    Event *event = new Event();
    Muon.Branch("events", &event);     
    TRandom2 *r2=new TRandom2();

    Initialize(1);
    for (int i = 0; i < Ngen; i++)
    {
        event->Reset();
        double moon_phi = 100. / 57.3;
        double moon_the = 35. / 57.3;
        double range_moon = 45. / 57.3 / 2.;
        double moon_size = 0.25 / 57.3;
        int flag = 0;
        do {
            flag = 0;
            the_gen = TMath::ACos(gRandom->Rndm());
            phi_gen = (2 * range_moon) * (gRandom->Rndm()) + (moon_phi - range_moon);
            double tmp = TMath::ACos(sin(moon_the) * sin(the_gen) * cos(phi_gen) * cos(moon_phi)
                    + sin(moon_the) * sin(the_gen) * sin(phi_gen) * sin(moon_phi)
                    + cos(moon_the) * cos(the_gen));
            if(tmp <= moon_size) flag = 1;
            else flag = 0;
        }while(flag);

        theta_phi_gen -> Fill(57.3 * (phi_gen - moon_phi) * TMath::Sin(moon_the), 57.3 * (the_gen - moon_the));
        phi->Fill(57.3*phi_gen);
        theta->Fill(57.3*the_gen);
        event->SetPhiGen(phi_gen);
        event->SetThetaGen(the_gen);
        
        cos_gen[0] = sin(the_gen)*cos(phi_gen);
        cos_gen[1] = sin(the_gen)*sin(phi_gen);
        cos_gen[2] = cos(the_gen);
        x_coord = 1200. * (gRandom -> Rndm()) - 600.;
        y_coord = 1200. * (gRandom -> Rndm()) - 600.;

        nstn = 0;
        for (int j = 0; j < N_STN_TOT; j++) {
            double R = sqrt((MuonX[j] - x_coord)*(MuonX[j] - x_coord) + (MuonY[j] - y_coord)* (MuonY[j] - y_coord));
            if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) {
                Station st;
                st.SetID(j + 1);
                st.SetX(MuonX[j]);
                st.SetY(MuonY[j]);
                st.SetZ(MuonZ[j]);
                st.SetT((cos_gen[0]*MuonX[j]  + cos_gen[1]*MuonY[j] + cos_gen[2]*MuonZ[j]) / 0.3 + TIME_ERROR * (r2->Gaus(0.,1.)));
                st.SetTerr(MuonTerr[j]);
                event->AddStation(st);
                nstn++;
            }
        }

        if(nstn > 2)
        {
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
        else i--;
    }
    auto c4 = new TCanvas("c4","");
    phi -> Draw();
    auto c5 = new TCanvas("c5","");
    theta -> Draw();
    auto c6 = new TCanvas("c6","");
    theta_phi_gen -> Draw("CONT4Z");
    return 0;
}
