#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;
#include <stdlib.h>
#include <bitset>
/*#if defined(__CINT__) && !defined(__MAKECINT__)
#include "AutoDict_Event_cxx.so"
#else
#include "Tunka.h"
#endif*/

static const int    N_STN_TOT =  19;
static const double P_STN_TRG = 1.;
static const double R0 = 150.; //meters
static const double COS_THE_POWER = 6.;
const char*         FILE_NAME    = "analysis_results/ToyMCRecon.root";
double              sep = 200000.;
int counterNoise = 0;
int counterSignal = 0;

vector<float> MuonX;
vector<float> MuonY;
vector<float> MuonZ;
vector<float> MuonTerr;


TH1D *cos_the_distrHISC = new TH1D("cos_the_distrHISC", "cos(theta)_HiSCORE", 500., 0., 1.1);
TH1D *cos_the_distrHISC_MC = new TH1D("cos_the_distrHISC_MC", "cos(theta)_HiSCORE", 500., 0., 1.1);

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
double FindMaxTimeInEvt(Event *evt)
{
    double maxTime = evt->MuStn(0).T();
    //cout << "Times: "<< evt->MuStn(0).T()<< endl;
    for(int i = 1; i < evt->NStations(); i++)
    {
       // cout << "Times: "<< evt->MuStn(i).T()<< endl;
        if(evt->MuStn(i).T() > maxTime) maxTime = evt->MuStn(i).T();
    }
   // cout << "MAX: "<< maxTime << "\n" << endl;
    return maxTime;
}
//===============================================================================================//
double FindMinTimeInEvt(Event *evt)
{
    double minTime = evt->MuStn(0).T();
    for(int i = 1; i < evt->NStations(); i++)
    {
        if(evt->MuStn(i).T() < minTime) minTime = evt->MuStn(i).T();
    }
    return minTime;
}
//===============================================================================================//
//Events generation and writing to the root-file
int Generate(int Ngen = 1000, const char* fileToSave = FILE_NAME)
{
    TH1D *phi = new TH1D("phi",   "phi", 370, -2,  362);
    TH1D *theta = new TH1D("theta",   "theta", 360,  -180,  180);
    TH2D *points = new TH2D("points", "phi sinHiscTheta, theta", 360, -180, 180, 90, -90, 90);

    int nstn;
    int syst = 0.;
    double the_gen, phi_gen, cos_gen[3], x_coord, y_coord;
    TTree Muon("Muon","");
    Event *event = new Event();
    Muon.Branch("events", &event);     
    TRandom2 *r2=new TRandom2();

    double Nsignal = 0.;

    Initialize(1);
    for (int i = 0; i < Ngen; i++)
    {
        event->Reset();
        int flag = 0;
        double HiSC_the = 20. / (180. / TMath::Pi());
        double HiSC_view_size = 30. / (180. / TMath::Pi());


        do {
            the_gen = TMath::ACos(pow(gRandom->Rndm(), 1. / (1. + COS_THE_POWER)));
            phi_gen = TMath::Pi() * (2 * gRandom->Rndm() - 1);

            if(TMath::ACos(sin(HiSC_the) * sin(the_gen) * cos(phi_gen) + cos(HiSC_the) * cos(the_gen))
            <= HiSC_view_size) flag = 0;
            else flag = 1;
        }while(flag);

        double a = phi_gen, b = the_gen;
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
        event->SetPhiGen(a);
        event->SetThetaGen(b);

        cos_gen[0] = sin(b)*cos(a);
        cos_gen[1] = sin(b)*sin(a);
        cos_gen[2] = cos(b);
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
                double timeGen = (cos_gen[0]*MuonX[j]  + cos_gen[1]*MuonY[j] + cos_gen[2]*MuonZ[j]) / 0.3 + MuonTerr[j] *(r2->Gaus(0.,1.));
                if(j + 1 == 1) st.SetT(timeGen + syst);
                else st.SetT(timeGen);

                int trow[4];
                trow[0] = 0;
                trow[1] = 0;
                trow[2] = 0;
                trow[3] = (int) timeGen;
                st.SetTrow(trow);

                st.SetTerr(MuonTerr[j]);
                event->AddStation(st);
                event->SetMatchHiSC(1);
                nstn++;
            }
        }

   /*  //Noise generation
        //if(i > Ngen * 0.46 && nstn >= 2)
            if(i > Ngen * 0.5 && nstn == 2)
            {
                int flag = 0;
                int id_noise;
                do{
                    flag = 0;
                    id_noise = gRandom -> Rndm() * (N_STN_TOT - 1);
                    for(int ii = 0; ii < nstn; ii++)
                    {
                        if(id_noise == event->MuStn(ii).ID()) flag = 1;
                    }
                }while(flag);

                //double t_noise = 2 * 10000. * gRandom -> Rndm() - 10000.;
                double t_noise = 10000. * gRandom -> Rndm();

                Station st;
                st.SetID(id_noise);
                st.SetX(MuonX[id_noise]);
                st.SetY(MuonY[id_noise]);
                st.SetZ(MuonZ[id_noise]);
                for(int ii = 0; ii < nstn; ii++) {
                    cout << event->MuStn(ii).T() << "\t" ;
                }
                cout << t_noise << endl;
                int trow[4];
                trow[0] = 0;
                trow[1] = 0;
                trow[2] = 0;
                trow[3] = (int) t_noise;
                st.SetTrow(trow);
                st.SetT(t_noise);

                st.SetTerr(MuonTerr[id_noise]);
                event->AddStation(st);
                event->SetMatchHiSC(1);
                nstn++;
                if(nstn == 3) counterNoise++;
            }
               //else if(i > Ngen && nstn != 2) i--;
       //end of noise generation*/


        if(nstn >= 3) Nsignal++;

        if(nstn >= 3)
        {
            if(nstn == 3) counterSignal++;
            points -> Fill(180. / TMath::Pi() * event->PhiGen() * TMath::Sin(HiSC_the), 180. / TMath::Pi() * event -> ThetaGen());
            Muon.Fill();
            cos_the_distrHISC_MC -> Fill(TMath::Cos(event->ThetaGen()));
            if(i == Ngen - 1)
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
    auto c9 = new TCanvas("c9","");
    points -> Draw("CONT4Z");
    cout << "\tNnoise part = " << counterNoise * 1. /  (counterSignal * 1. - counterNoise * 1.) * 100. << endl;
    return 0;
}
//======================================================================================================================
//For drawing MC and data angular distr
void FillHists(const char *FileName, int draw = 1)
{
    vector<Event> evts;
    ReadTreeOfEvts(FileName, evts);
    for(int i = 0; i < (int) evts.size(); i++)
    {
        if(evts[i].MatchHiSC() == 1 && evts[i].NStations() > 2 && evts[i].Chi2() < 10)
        {
            cos_the_distrHISC -> Fill(TMath::Cos(evts[i].ThetaGen()));
        }
    }
    if(draw)
    {
        auto c11 = new TCanvas("c11", "");
        cos_the_distrHISC -> Scale(1. / cos_the_distrHISC -> GetEntries());
        cos_the_distrHISC_MC -> Scale(1. / cos_the_distrHISC_MC -> GetEntries());
        cos_the_distrHISC -> GetXaxis() -> SetTitle("cos(theta)");
        cos_the_distrHISC -> Draw();
        cos_the_distrHISC_MC -> SetLineColor(kRed);
        cos_the_distrHISC_MC -> Draw("same hist");
    }

}
 int DrawCompare()
 {
    Generate(100000, "toyMC.root");

    FillHists("analysis_results/for_EAS_center/Gr+Hisc221119Recon.root", 0);
    FillHists("analysis_results/for_EAS_center/Gr+Hisc261119Recon.root", 0);
    FillHists("analysis_results/for_EAS_center/Gr+Hisc271119Recon.root", 0);
    FillHists("analysis_results/for_EAS_center/Gr+Hisc301219Recon.root", 0);
    FillHists("analysis_results/for_EAS_center/Gr+Hisc150220Recon.root", 0);
    FillHists("analysis_results/for_EAS_center/Gr+Hisc160220Recon.root", 0);
    FillHists("analysis_results/for_EAS_center/Gr+Hisc200220Recon.root", 1);
    return 0;
 }