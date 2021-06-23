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
static const double P_STN_TRG = 1.;
static const double R0 = 150.; //meters
static const double COS_THE_POWER = 6.;
const char*         FILE_NAME    = "analysis_results/ToyMC.root";

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
int Generate(int Ngen = 1000)
{
    TH1D *phi = new TH1D("phi",   "phi", 370, -2,  362);
    TH1D *theta = new TH1D("theta",   "theta", 360,  -180,  180);
    TH1D *angleBtwEvts = new TH1D("angleBtwEvts", "angle between 2st events", 100., -3, 181.);

    int nstn;
    double the_gen, phi_gen, cos_gen[3], x_coord, y_coord, the_gen2, phi_gen2;
    TRandom2 *r2=new TRandom2();

    Initialize(1);
    for (int i = 0; i < Ngen; i++)
    {
       // int flag = 0;
        double HiSC_the = 20. / (180. / TMath::Pi());
        double HiSC_view_size = 30. / (180. / TMath::Pi());

     //   do {
            the_gen = TMath::ACos(pow(gRandom->Rndm(), 1. / 7.));
            phi_gen = TMath::Pi() * (2 * gRandom->Rndm() - 1);

            the_gen2 = TMath::ACos(pow(gRandom->Rndm(), 1. / 7.));
            phi_gen2 = TMath::Pi() * (2 * gRandom->Rndm() - 1);

            //double tmp = the_gen - HiSC_the;
          //  double tmp2 = phi_gen;

           // if(TMath::ACos(sin(HiSC_the) * sin(the_gen) * cos(phi_gen) + cos(HiSC_the) * cos(the_gen))
           // <= HiSC_view_size) flag = 0;
          //  else flag = 1;
    //    }while(flag);

        double a = phi_gen, b = the_gen;
        double a2 = phi_gen2, b2 = the_gen2;
        if(b < 0)
        {
            b = -b;
            if(a > 0) a -= TMath::Pi();
            else a += TMath::Pi();
        }
        if(a < 2 * TMath::Pi()) a += 2 * TMath::Pi();
        if(a >= 2 * TMath::Pi()) a -= 2 * TMath::Pi();
        if(a < 0) a += 2 * TMath::Pi();
        phi->Fill(180. / TMath::Pi() * a);
        theta->Fill(180. / TMath::Pi() * b);

        cos_gen[0] = sin(b)*cos(a);
        cos_gen[1] = sin(b)*sin(a);
        cos_gen[2] = cos(b);
        x_coord = 1200. * (gRandom -> Rndm()) - 600.;
        y_coord = 1200. * (gRandom -> Rndm()) - 600.;
        nstn = 0;
        for (int j = 0; j < N_STN_TOT; j++) {
            double R = sqrt((MuonX[j] - x_coord)*(MuonX[j] - x_coord) + (MuonY[j] - y_coord)* (MuonY[j] - y_coord));
            if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) nstn++;
        }

        if(nstn == 2)
        {
            double dpsi = 57.3 * TMath::ACos(cos(b) * cos(b2) + sin(b) * sin(b2) * cos(a) * cos(a2) +
                                             sin(b) * sin(b2) * sin(a) * sin(a2));
            angleBtwEvts -> Fill(dpsi);
        }
        else i--;
    }

    auto c4 = new TCanvas("c4","");
    phi -> Draw();
    auto c5 = new TCanvas("c5","");
    theta -> Draw();
    auto c9 = new TCanvas("c9","");
    angleBtwEvts -> Draw();
    return 0;
}
