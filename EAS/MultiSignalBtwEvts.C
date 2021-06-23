//Find multisignal events
#include "ReadTreeOfEvts.h"
using namespace std;

static const int TIME_WINDOW = 125; // s
static const int N_STN_MIN =   3;
const char      *FILE_NAME = "analysis_results/Gr+Hisc241119Recon.root";
TH1D *angleBtwEvtsExp2 = new TH1D("angleBtwEvtsExp2", "dpsi_exp 2 evts", 80., -3, 181.);
TH1D *angleBtwEvtsExp3 = new TH1D("angleBtwEvtsExp3", "dpsi_exp 3 evts", 60., -3, 181.);
TH1D *angleBtwEvtsMC2 = new TH1D("angleBtwEvtsMC2", "dpsi_toyMC 2 evts", 80., -3, 181.);
TH1D *angleBtwEvtsMC3 = new TH1D("angleBtwEvtsMC3", "dpsi_toyMC 3 evts", 60., -3, 181.);
//TH1D *signum_ = new TH1D("signum", "", 4000., -2000000, 2000000.);
TH2F *timeDiff = new TH2F("timeDiff", "time difference", 210., -10, 200., 210., -10, 200.);

static const int    N_STN_TOT =  19;
static const double P_STN_TRG = 1.;
static const double R0 = 150.; //meters
static const double COS_THE_POWER = 6.;

vector<float> MuonX;
vector<float> MuonY;
vector<float> MuonZ;
vector<float> MuonTerr;
//----------------------------------------------------------------------------------------------------------------------
void CountNstEvents(vector<Event> &evts)
{
    int counter3 = 0, counter4 = 0;
    for(int i = 0; i < (int)evts.size(); i++)
    {
        if(evts[i].NStations()  >= 3) counter3++;
        if(evts[i].NStations() >= 4) counter4++;
    }
    cout << counter3 << " - 3st evts\t" << counter4 <<" - 4st evts" << endl;
}
//----------------------------------------------------------------------------------------------------------------------
void SelectNstMinEvts(vector<Event> &evts0, vector<Event> &evts)
{
    for(int i = 0; i < (int)evts0.size(); i++)
    {
        if(evts0[i].NStations() >= N_STN_MIN) evts.push_back(evts0[i]);
    }
    printf("3 st evts are selected\n");
}
//----------------------------------------------------------------------------------------------------------------------
float GetTimeNs(Event &evt)
{
    float t = evt.EvtTime()[0] * 3600. + evt.EvtTime()[1] * 60. + evt.EvtTime()[2];
    t = t * 1000000000. + evt.EvtTime()[3];
    return t;
}
//----------------------------------------------------------------------------------------------------------------------
void ReculcEvtsTime(vector<Event> &evts)
{
    for(int i = 0; i < (int)evts.size(); i++)
    {
        double t0 = evts[i].MuStn(0).Trow()[0] * 3600. + evts[i].MuStn(0).Trow()[1] * 60. + evts[i].MuStn(0).Trow()[2];
        t0 = t0 * 1000000000. + evts[i].MuStn(0).Trow()[3];
        int j_min = 0;
        for(int j = 1; j < evts[i].NStations(); j++)
        {
            double t = evts[i].MuStn(j).Trow()[0] * 3600. + evts[i].MuStn(j).Trow()[1] * 60. + evts[i].MuStn(j).Trow()[2];
            t = t * 1000000000. + evts[i].MuStn(j).Trow()[3];
            if(t < t0){
                t0 = t;
                j_min = j;
            }
        }
        evts[i].SetEvtTime(evts[i].MuStn(j_min).Trow());
    }
}
//----------------------------------------------------------------------------------------------------------------------
float DtMin(vector<Event> &newVect, double t2)
{
    float dt = t2 - GetTimeNs(newVect[0]);
    for (int l = 1; l < (int) newVect.size(); l++)
    {
        float t1 = GetTimeNs(newVect[l]);
        float dt_ = t2 - t1;
        if (fabs(dt_) < fabs(dt)) dt = dt_;
    }
    return dt;
}
//----------------------------------------------------------------------------------------------------------------------
float MaxT(vector<Event> &newVect)
{
    float t = GetTimeNs(newVect[0]);
    for(int l = 0; l < (int) newVect.size(); l++)
    {
        float t0 = GetTimeNs(newVect[l]);
        if (t0 > t) t = t0;
    }
    return t;
}
//----------------------------------------------------------------------------------------------------------------------
/*void FillNewEvtsVector(vector<Event> &evts, vector<vector<Event>> &newEvts)
{
    printf("start to select multiple events\n");
    for(int i = 0; i < (int)evts.size(); i++)
    {
        vector<Event> newVect;
        newVect.push_back(evts[i]);
        for(int j = i + 1; j < (int)evts.size(); j++)
        {
            double time2 = GetTimeNs(evts[j]);
            double time1 = MaxT(newVect);
            double dt = time2 - time1;
            //double dt = DtMin(newVect, time2);
            if (fabs(dt) <= TIME_WINDOW) newVect.push_back(evts[j]);
            else if (dt > 0 && fabs(dt) > TIME_WINDOW) break;
        }
        if((int)newVect.size() > 1)
        {
            newEvts.push_back(newVect);
            newVect.clear();
        }
    }
    cout << newEvts.size()  << " multi-events are found from " << evts.size() << " evts\n" << endl;
}*/
void FillNewEvtsVector(vector<Event> &evts, vector<vector<Event>> &newEvts)
{
    printf("start to select multiple events\n");
    for(int i = 0; i < (int)evts.size(); i++)
    {
        vector<Event> newVect;
        newVect.push_back(evts[i]);
        float time1 = GetTimeNs(evts[i]);
        for(int j = 0; j < (int)evts.size(); j++)
        {
            float time2 = GetTimeNs(evts[j]);
            float dt = time2 - time1;
           // if(j - i < 4) signum_-> Fill(dt);
            if (dt > 0 && fabs(dt) < TIME_WINDOW * 1000000000.) newVect.push_back(evts[j]);
            else if (dt > 0 && dt > TIME_WINDOW * 1000000000.) break;
        }
        if((int)newVect.size() > 1)
        {
            newEvts.push_back(newVect);
            newVect.clear();
        }
    }
    cout << newEvts.size()  << " multi-events are found from " << evts.size() << " evts\n" << endl;
}
//----------------------------------------------------------------------------------------------------------------------
void FillHists2evt(vector<vector<Event>> &newEvts)
{
    printf("filling 2 evt histograms\n");
    for (int i = 0; i < (int) newEvts.size(); i++)
    {
        if ((int) newEvts[i].size() == 2)
        {
            double the0, the1, phi0, phi1;
            the0 = newEvts[i][0].ThetaRec();
            the1 = newEvts[i][1].ThetaRec();
            phi0 = newEvts[i][0].PhiRec();
            phi1 = newEvts[i][1].PhiRec();
            double dpsi = 57.3 * TMath::ACos(cos(the0) * cos(the1) + sin(the0) * sin(the1) * cos(phi0) * cos(phi1) +
                                             sin(the0) * sin(the1) * sin(phi0) * sin(phi1));
            angleBtwEvtsExp2->Fill(dpsi);

            if (newEvts[i][0].ThetaGen() != -100) {
                the0 = newEvts[i][0].ThetaGen();
                the1 = newEvts[i][1].ThetaGen();
                phi0 = newEvts[i][0].PhiGen();
                phi1 = newEvts[i][1].PhiGen();
                double dpsi = 57.3 * TMath::ACos(cos(the0) * cos(the1) + sin(the0) * sin(the1) * cos(phi0) * cos(phi1) +
                                                 sin(the0) * sin(the1) * sin(phi0) * sin(phi1));
                //cout << i << " \t "<< newEvts[i].size() << " evts; dpsiHiSC = " << dpsi << endl;
                float dx = fabs(newEvts[i][0].EAScentersGen()[0] - newEvts[i][1].EAScentersGen()[0]);
                float dy = fabs(newEvts[i][0].EAScentersGen()[1] - newEvts[i][1].EAScentersGen()[1]);
                //cout << "dx, dy center HiSCORE: " << dx << "\t" << dy << "\n" << endl;
            }
        }

    }
}
//----------------------------------------------------------------------------------------------------------------------
double CalcAngleBtwDirections(double phi1, double phi2, double the1, double the2)
{
    return (180. / TMath::Pi() * TMath::ACos(cos(the2) * cos(the1) + sin(the2) * sin(the1) * cos(phi2) * cos(phi1) +
                              sin(the2) * sin(the1) * sin(phi2) * sin(phi1)));
}
//----------------------------------------------------------------------------------------------------------------------
void FillHists3evt(vector<vector<Event>> &newEvts)
{
    printf("filling 3 evt histograms\n");
    for(int i = 0; i < (int)newEvts.size(); i++)
    {
        if((int)newEvts[i].size() == 3)
        {
            double the0, the1, phi0, phi1, the2, phi2;
            the0 = newEvts[i][0].ThetaRec();
            the1 = newEvts[i][1].ThetaRec();
            the2 = newEvts[i][2].ThetaRec();
            phi0 = newEvts[i][0].PhiRec();
            phi1 = newEvts[i][1].PhiRec();
            phi2 = newEvts[i][2].PhiRec();
            double dpsi0 = CalcAngleBtwDirections(phi0, phi1, the0, the1);
            double dpsi1 = CalcAngleBtwDirections(phi0, phi2, the0, the2);
            double dpsi2 = CalcAngleBtwDirections(phi1, phi2, the1, the2);
            angleBtwEvtsExp3 -> Fill((dpsi0 + dpsi1 + dpsi2) / 3.);

            float dt[2], t[3];
            for(int u = 0; u < 3; u++) t[u] = GetTimeNs(newEvts[i][u]);
            dt[0] = fabs((t[1] - t[0])) / 1000000000.;
            dt[1] = fabs((t[2] - t[0])) / 1000000000.;
            //dt[2] = fabs((t[2] - t[1])) / 1000000000.;
            sort(dt, dt + 2);
            timeDiff -> Fill(dt[0] , dt[1]);

           // if(dt[0] > TIME_WINDOW / 1000000000. || dt[1] > TIME_WINDOW / 1000000000.  || dt[2] > TIME_WINDOW / 1000000000. )cout<<i << "\t" << t[0] << ":"
            //<< t[1] << ":" << t[2] << "\n" << dt[0] << "\t" <<dt[1] << "\t" << dt[2] << "\n\n";
        }
    }
}
//----------------------------------------------------------------------------------------------------------------------
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
int Generate2(int Ngen = 1000, int poinSource = 0)
{
    int nstn1, nstn2, nstn3;
    double the_gen, phi_gen, cos_gen[3], x_coord, y_coord, the_gen2, phi_gen2;
    TRandom2 *r2=new TRandom2();

    for (int i = 0; i < Ngen; i++)
    {
        the_gen = TMath::ACos(pow(gRandom->Rndm(), 1. / 7.));
        phi_gen = TMath::Pi() * (2 * gRandom->Rndm() - 1);
        the_gen2 = TMath::ACos(pow(gRandom->Rndm(), 1. / 7.));
        phi_gen2 = TMath::Pi() * (2 * gRandom->Rndm() - 1);

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

        if(b2 < 0)
        {
            b2 = -b2;
            if(a2 > 0) a2 -= TMath::Pi();
            else a2 += TMath::Pi();
        }
        if(a2 < 2 * TMath::Pi()) a2 += 2 * TMath::Pi();
        if(a2 >= 2 * TMath::Pi()) a2 -= 2 * TMath::Pi();
        if(a2 < 0) a2 += 2 * TMath::Pi();

        cos_gen[0] = sin(b)*cos(a);
        cos_gen[1] = sin(b)*sin(a);
        cos_gen[2] = cos(b);
        x_coord = 1200. * (gRandom -> Rndm()) - 600.;
        y_coord = 1200. * (gRandom -> Rndm()) - 600.;
        nstn1 = 0;
        for (int j = 0; j < N_STN_TOT; j++) {
            double R = sqrt((MuonX[j] - x_coord)*(MuonX[j] - x_coord) + (MuonY[j] - y_coord)* (MuonY[j] - y_coord));
            if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) nstn1++;
        }
        cos_gen[0] = sin(b2)*cos(a2);
        cos_gen[1] = sin(b2)*sin(a2);
        cos_gen[2] = cos(b2);
        x_coord = 1200. * (gRandom -> Rndm()) - 600.;
        y_coord = 1200. * (gRandom -> Rndm()) - 600.;
        nstn2 = 0;
        for (int j = 0; j < N_STN_TOT; j++) {
            double R = sqrt((MuonX[j] - x_coord)*(MuonX[j] - x_coord) + (MuonY[j] - y_coord)* (MuonY[j] - y_coord));
            if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) nstn2++;
        }

        if(nstn1 >= N_STN_MIN && nstn2 >= N_STN_MIN)
        {
            double dpsi = 57.3 * TMath::ACos(cos(b) * cos(b2) + sin(b) * sin(b2) * cos(a) * cos(a2) +
                                             sin(b) * sin(b2) * sin(a) * sin(a2));
            angleBtwEvtsMC2 -> Fill(dpsi);
        }
        else i--;
    }

    if(poinSource == 1)
    {
        for(int i = 0; i < 1000; i++)
        {
            the_gen = 70. / 57.3;
            phi_gen = 0.;
            the_gen2 = 70. / 57.3;
            phi_gen2 = 0.;
            double a = phi_gen, b = the_gen;
            double a2 = phi_gen2, b2 = the_gen2;

            cos_gen[0] = sin(b)*cos(a);
            cos_gen[1] = sin(b)*sin(a);
            cos_gen[2] = cos(b);
            x_coord = 1200. * (gRandom -> Rndm()) - 600.;
            y_coord = 1200. * (gRandom -> Rndm()) - 600.;
            nstn1 = 0;
            for (int j = 0; j < N_STN_TOT; j++) {
                double R = sqrt((MuonX[j] - x_coord)*(MuonX[j] - x_coord) + (MuonY[j] - y_coord)* (MuonY[j] - y_coord));
                if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) nstn1++;
            }
            cos_gen[0] = sin(b2)*cos(a2);
            cos_gen[1] = sin(b2)*sin(a2);
            cos_gen[2] = cos(b2);
            x_coord = 1200. * (gRandom -> Rndm()) - 600.;
            y_coord = 1200. * (gRandom -> Rndm()) - 600.;
            nstn2 = 0;
            for (int j = 0; j < N_STN_TOT; j++) {
                double R = sqrt((MuonX[j] - x_coord)*(MuonX[j] - x_coord) + (MuonY[j] - y_coord)* (MuonY[j] - y_coord));
                if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) nstn2++;
            }

            if(nstn1 >= N_STN_MIN && nstn2 >= N_STN_MIN) angleBtwEvtsMC2 -> Fill(0.);
            else i--;
        }
    }
    return 0;
}
//===============================================================================================//
//Events generation and writing to the root-file
int Generate3(int Ngen = 1000)
{
    int nstn1, nstn2, nstn3;
    double the_gen, phi_gen, cos_gen[3], x_coord, y_coord, the_gen2, phi_gen2, the_gen3, phi_gen3;
    TRandom2 *r2=new TRandom2();

    for (int i = 0; i < Ngen; i++)
    {
        the_gen = TMath::ACos(pow(gRandom->Rndm(), 1. / 7.));
        phi_gen = TMath::Pi() * (2 * gRandom->Rndm() - 1);
        the_gen2 = TMath::ACos(pow(gRandom->Rndm(), 1. / 7.));
        phi_gen2 = TMath::Pi() * (2 * gRandom->Rndm() - 1);
        the_gen3 = TMath::ACos(pow(gRandom->Rndm(), 1. / 7.));
        phi_gen3 = TMath::Pi() * (2 * gRandom->Rndm() - 1);

        double a = phi_gen, b = the_gen;
        double a2 = phi_gen2, b2 = the_gen2;
        double a3 = phi_gen3, b3 = the_gen3;
        if(b < 0)
        {
            b = -b;
            if(a > 0) a -= TMath::Pi();
            else a += TMath::Pi();
        }
        if(a < 2 * TMath::Pi()) a += 2 * TMath::Pi();
        if(a >= 2 * TMath::Pi()) a -= 2 * TMath::Pi();
        if(a < 0) a += 2 * TMath::Pi();

        if(b2 < 0)
        {
            b2 = -b2;
            if(a2 > 0) a2 -= TMath::Pi();
            else a2 += TMath::Pi();
        }
        if(a2 < 2 * TMath::Pi()) a2 += 2 * TMath::Pi();
        if(a2 >= 2 * TMath::Pi()) a2 -= 2 * TMath::Pi();
        if(a2 < 0) a2 += 2 * TMath::Pi();

        if(b3 < 0)
        {
            b3 = -b3;
            if(a3 > 0) a3 -= TMath::Pi();
            else a3 += TMath::Pi();
        }
        if(a3 < 2 * TMath::Pi()) a3 += 2 * TMath::Pi();
        if(a3 >= 2 * TMath::Pi()) a3 -= 2 * TMath::Pi();
        if(a3 < 0) a3 += 2 * TMath::Pi();

        cos_gen[0] = sin(b)*cos(a);
        cos_gen[1] = sin(b)*sin(a);
        cos_gen[2] = cos(b);
        x_coord = 1200. * (gRandom -> Rndm()) - 600.;
        y_coord = 1200. * (gRandom -> Rndm()) - 600.;
        nstn1 = 0;
        for (int j = 0; j < N_STN_TOT; j++) {
            double R = sqrt((MuonX[j] - x_coord)*(MuonX[j] - x_coord) + (MuonY[j] - y_coord)* (MuonY[j] - y_coord));
            if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) nstn1++;
        }
        cos_gen[0] = sin(b2)*cos(a2);
        cos_gen[1] = sin(b2)*sin(a2);
        cos_gen[2] = cos(b2);
        x_coord = 1200. * (gRandom -> Rndm()) - 600.;
        y_coord = 1200. * (gRandom -> Rndm()) - 600.;
        nstn2 = 0;
        for (int j = 0; j < N_STN_TOT; j++) {
            double R = sqrt((MuonX[j] - x_coord)*(MuonX[j] - x_coord) + (MuonY[j] - y_coord)* (MuonY[j] - y_coord));
            if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) nstn2++;
        }
        cos_gen[0] = sin(b3)*cos(a3);
        cos_gen[1] = sin(b3)*sin(a3);
        cos_gen[2] = cos(b3);
        x_coord = 1200. * (gRandom -> Rndm()) - 600.;
        y_coord = 1200. * (gRandom -> Rndm()) - 600.;
        nstn3 = 0;
        for (int j = 0; j < N_STN_TOT; j++) {
            double R = sqrt((MuonX[j] - x_coord)*(MuonX[j] - x_coord) + (MuonY[j] - y_coord)* (MuonY[j] - y_coord));
            if ((gRandom -> Rndm()) < P_STN_TRG * exp(-R/R0)) nstn3++;
        }

        if(nstn1 >= N_STN_MIN && nstn2 >= N_STN_MIN && nstn3 >+ N_STN_MIN)
        {
            double dpsi1 = 57.3 * TMath::ACos(cos(b) * cos(b2) + sin(b) * sin(b2) * cos(a) * cos(a2) +
                                              sin(b) * sin(b2) * sin(a) * sin(a2));
            double dpsi2 = 57.3 * TMath::ACos(cos(b) * cos(b3) + sin(b) * sin(b3) * cos(a) * cos(a3) +
                                              sin(b) * sin(b3) * sin(a) * sin(a3));
            double dpsi3 = 57.3 * TMath::ACos(cos(b3) * cos(b2) + sin(b3) * sin(b2) * cos(a3) * cos(a2) +
                                              sin(b3) * sin(b2) * sin(a3) * sin(a2));
            angleBtwEvtsMC3 -> Fill((dpsi1 + dpsi2 + dpsi3) / 3.);
        }
        else i--;
    }
    return 0;
}
//----------------------------------------------------------------------------------------------------------------------
void DrawHists()
{
    auto c100 = new TCanvas("c100","");
    timeDiff -> SetMarkerStyle(7);
    timeDiff -> Draw();

    Initialize(1);

    printf("simulation 2 evts\n");
    Generate2(20000, 0);
    printf("simulation 3 evts\n");
    Generate3(20000);

    auto c9 = new TCanvas("c9","");
    angleBtwEvtsExp2 -> Scale(1. / angleBtwEvtsExp2 -> GetEntries());
    angleBtwEvtsExp2 -> Draw();
    angleBtwEvtsMC2 -> Scale(1. / angleBtwEvtsMC2 -> GetEntries());
    angleBtwEvtsMC2 -> SetLineColor(kRed);
    angleBtwEvtsMC2 -> Draw("SAME HIST");
    auto c15 = new TCanvas("c15","");
    angleBtwEvtsExp3 -> Scale(1. / angleBtwEvtsExp3 -> GetEntries());
    angleBtwEvtsExp3 -> Draw();
    angleBtwEvtsMC3 -> SetLineColor(kRed);
    angleBtwEvtsMC3 -> Scale(1. / angleBtwEvtsMC3 -> GetEntries());
    angleBtwEvtsMC3 -> Draw("SAME HIST");

    printf("clear memory\n\n");
    MuonX.clear();
    MuonY.clear();
    MuonZ.clear();
    MuonTerr.clear();
}
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
int FindMultiSignalBtwEvts(const char * InputFileName = FILE_NAME, int draw = 1)
{
    vector<Event> evts0, evts;
    vector<vector<Event>> newEvts;
    printf("reading the file\n");
    ReadTreeOfEvts(InputFileName, evts0);

    CountNstEvents(evts0);
    SelectNstMinEvts(evts0, evts);
    evts0.clear();
    ReculcEvtsTime(evts);
    FillNewEvtsVector(evts, newEvts);//multiple evt vector
    FillHists2evt(newEvts);
    FillHists3evt(newEvts);
    if(draw)
    {
        DrawHists();
        evts.clear();
        newEvts.clear();
    }
    return 0;
}
//----------------------------------------------------------------------------------------------------------------------
int main3()
{
    FindMultiSignalBtwEvts("analysis_results/Gr+Hisc241119Recon.root", 0);
    FindMultiSignalBtwEvts("analysis_results/Gr+Hisc241219Recon.root", 0);
    FindMultiSignalBtwEvts("analysis_results/Gr+Hisc291219Recon.root", 0);
    FindMultiSignalBtwEvts("analysis_results/Gr+Hisc020120Recon.root", 0);
    FindMultiSignalBtwEvts("analysis_results/Gr+Hisc210220Recon.root");
   // auto c19 = new TCanvas("c19","");
   // signum_->Draw();
    return 0;
}