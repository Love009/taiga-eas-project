//Find multisignal inside one event time interval
#include "ReadTreeOfEvts.h"
using namespace std;

const int TIME_WINDOW = 5; // *5 ns
const int LOW_SIG_CHECK_INTERVAL = 30; // counts of 200 total ADS counts
static const int    N_STN_MIN =   2;
vector<double> MeanCH0;
vector<double> MeanCH2;
vector<double> SigmaCH0;
vector<double> SigmaCH2;
TH1D *signalStart = new TH1D("signalStart","signalStart", 200,  0,  200); //hist for linear fit
const char *INITIALIZATION_FILE_NAME = "analysis_results/pedestals/pedestal241219.txt";
const char *INPUT_FILE_NAME          = "analysis_results/Gr+Hisc241219.root";
vector<Event> evts;

//========================================================================================================================================//
//========================================================================================================================================//
//========================================================================================================================================//
//Initualize station channel pedestals
void Initialize(const char * FileName = INITIALIZATION_FILE_NAME)
{
    double meanCH0, meanCH2, sigmaCH0, sigmaCH2;
    string  title, rest;
    int i = 1;
    ifstream in(FileName);
    if(in.fail()){
        cout << "File " << FileName << " does not exist" << endl;
    }
    else {
        getline(in, title);
        while (1) {
            in >> rest >> meanCH0 >> sigmaCH0 >> meanCH2 >> sigmaCH2;
            if (!in.good()) break;
            MeanCH0.push_back(meanCH0);
            MeanCH2.push_back(meanCH2);
            SigmaCH0.push_back(sigmaCH0);
            SigmaCH2.push_back(sigmaCH2);
            i++;
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Find a signal start using linear fit (the least square method)
int FindSignal(vector<int> ADCsig, int stID, int chanID, double sigma, int threshold, int lowLimit, int window = TIME_WINDOW)
{
    double sigTime = -1., flag = 0;
    double mean = MeanCH0[stID - 1];
    if(chanID == 2) mean = MeanCH2[stID - 1];
    for(int k = 0; k < (int) ADCsig.size(); k++)
    {
        int delta = ADCsig[k] -  mean;
        if(delta <= lowLimit * sigma && k <= LOW_SIG_CHECK_INTERVAL)
        {
            cout << "Signal ampl is lower than " << lowLimit * sigma << "\tStation_number = " << stID << "\t channel " << chanID << "\n\n";
            // TCanvas *c1 = new TCanvas("c1", "Signal from channel");
            TH1I *sig = new TH1I("sig","Signal from channel 0, too low", (int)ADCsig.size(), 0, (int)ADCsig.size());
            for(int k = 0; k < (int)ADCsig.size(); k++) { sig->SetBinContent((k + 1), ADCsig[k]-MeanCH0[stID-1]);}
            //sig -> Draw();
            sigTime = -2.;
            break;
        }
        else if(flag >= window)
        {
            TArrayD x(2*window), y(window);
            TMatrixD A(window, 2), A_t(2, window),Y(window, 1), Pars(2, 1);
            for(int i = 0; i < window; i++)
            {
                y[i] = ADCsig[k - 1 - window + i] - mean;
                x[2*i] = k - 1 - window + i;
                x[2*i+1] = 1;
            }
            A.SetMatrixArray(x.GetArray());
            Y.SetMatrixArray(y.GetArray());
            A_t.Transpose(A);
            TMatrixD M = A_t*A;
            M.Invert();
            Pars = M * A_t * Y;
            double *pars = Pars.GetMatrixArray();
            sigTime = -pars[1] / pars[0];
            break;
        }
        else if(delta >= threshold * sigma) flag++;
        else flag = 0;
    }
    return (sigTime * 5);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CountVertices(vector<int> ADCsig, int stID, int chanID, double sigma, int threshold, int lowLimit, int window = TIME_WINDOW)
{
    int VerticesNum = 0, flag = 0, kMin = 0, sign = 1;
    double mean = MeanCH0[stID - 1];
    if(chanID == 2) mean = MeanCH2[stID - 1];
    for(int k = 0; k < (int) ADCsig.size(); k++)
    {
        if(VerticesNum == 0)
        {
            double delta = ADCsig[k] - mean;
            if (delta <= lowLimit * sigma && k <= LOW_SIG_CHECK_INTERVAL)
            {
               // cout << "Signal ampl is lower than " << lowLimit * sigma << "\tStation_number = " << stID
               //      << "\t channel " << chanID << "\n\n";
                break;
            }
            else if (flag >= window)
            {
                VerticesNum++;
                flag = 0;
                //break;
            }
            else if (delta >= threshold * sigma) flag++;
            else flag = 0;
        }
        else
        {
            double delta = 0.;
            for(int l = 0; l < 5; l++) delta += (ADCsig[k + l] - ADCsig[k - 1]) / (l + 1.);
            delta = delta / 5.;
            if(abs(delta) >= threshold * sigma)
            {
                if(delta / abs(delta) == sign) k += 4;
                else if(delta < 0 && sign > 0) {k += 4; sign = -1;}
                else if (delta > 0 && sign < 0 && kMin == 0) { kMin = k; sign = 1; k += 4;}
                else if (delta > 0 && kMin != 0) { VerticesNum++; kMin = 0; sign = 1; k += 4;}
            }
        }
    }
    return VerticesNum;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int DrawSepEvt(vector<Event> evts, int evtNum, int stNum) {
    cout << evts[evtNum].NStations() << " stations in event" << endl;

    TCanvas *c3 = new TCanvas("c3", "Signal from ADC channel");
    c3->Divide(1, 2);
    TH1I *sig0 = new TH1I("sig0", "Signal from ADC channel 0", (int) evts[evtNum].MuStn(stNum).ADC0().size(), 0,
                          (int) evts[evtNum].MuStn(stNum).ADC0().size());
    int stID = evts[evtNum].MuStn(stNum).ID();
    for (int k = 0; k < (int) evts[evtNum].MuStn(stNum).ADC0().size(); k++) { sig0->SetBinContent((k + 1), evts[evtNum].MuStn(
                stNum).ADC0()[k] - MeanCH0[stID - 1]);
    }
    c3->cd(1);
    sig0->Draw("");

    TH1I *sig2 = new TH1I("sig2", "Signal from ADC channel 2", (int) evts[evtNum].MuStn(stNum).ADC2().size(), 0,
                          (int) evts[evtNum].MuStn(stNum).ADC2().size());
    for (int k = 0; k < (int) evts[evtNum].MuStn(stNum).ADC2().size(); k++) { sig2->SetBinContent((k + 1), evts[evtNum].MuStn(
                stNum).ADC2()[k] - MeanCH2[stID - 1]);
    }
    c3->cd(2);
    sig2->Draw("");
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D *hist = new TH1D("hist", "time difference for signals of one station", 1000., -10., 2000.);
double FindDeadTime(vector<Event> evts, int stID)
{
    double deadTime = 0., tmp = 0.;
    for(int i = 0; i < (int)evts.size(); i++)
    {
        if(evts[i].NStations() >= 3) {
            for (int j = 0; j < evts[i].NStations(); j++) {
                if (evts[i].MuStn(j).ID() == stID) {
                    if (tmp == 0) {
                        tmp = evts[i].MuStn(j).Trow()[0] * 3600 + evts[i].MuStn(j).Trow()[1] * 60 +
                              evts[i].MuStn(j).Trow()[2];
                        tmp = tmp + evts[i].MuStn(j).Trow()[3] / 1000000000.;
                    } else {
                        double t = evts[i].MuStn(j).Trow()[0] * 3600 + evts[i].MuStn(j).Trow()[1] * 60 +
                                   evts[i].MuStn(j).Trow()[2];
                        t = t + evts[i].MuStn(j).Trow()[3] / 1000000000.;
                        hist->Fill(t - tmp);
                        tmp = t;
                    }
                }

            }
        }
    }
    TCanvas *c8 = new TCanvas("c8", "");
    hist -> Draw();
    return deadTime;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FindMultiSignal(const char * InitFileName = INITIALIZATION_FILE_NAME, const char * FileName = INPUT_FILE_NAME, int threshold = 5, int lowLimit = -10)
{
//    vector<Event> evts;
    ReadTreeOfEvts(FileName, evts);
    Initialize(InitFileName);
    int sigTime0/*ns*/, sigTime2/*ns*/, NumUsefulEvts = 0;
    for(int i = 0; i < (int) evts.size(); i++)
    {
        if(evts[i].MatchHiSC() == 1)
        {
            NumUsefulEvts++;
            int NstsWithMultiSig = 0;
            for (int j = 0; j < evts[i].NStations(); j++)
            {
                int stID = evts[i].MuStn(j).ID();

                sigTime0 = FindSignal(evts[i].MuStn(j).ADC0(), stID, 0, SigmaCH0[stID - 1], threshold, lowLimit);
                sigTime2 = FindSignal(evts[i].MuStn(j).ADC2(), stID, 2, SigmaCH2[stID - 1], threshold, lowLimit);

                if (sigTime0 > 0 && sigTime2 > 0) {
                    int nSignals0 = CountVertices(evts[i].MuStn(j).ADC0(), stID, 0, SigmaCH0[stID - 1], threshold - 2,
                                                  lowLimit);
                    int nSignals2 = CountVertices(evts[i].MuStn(j).ADC2(), stID, 2, SigmaCH2[stID - 1], threshold - 2,
                                                  lowLimit);
                    if ((nSignals0 > 1 || nSignals2 > 1) && NstsWithMultiSig != 0) cout << "\n\nStation num = " << j << "\n";
                    if (nSignals0 > 1 || nSignals2 > 1) NstsWithMultiSig++;
                }
            }
            if (NstsWithMultiSig >= 2 && evts[i].NStations() == 3) cout << "I've found multisignal (" << NstsWithMultiSig << " signals) in event  number " << i << endl;
        }
    }
    cout << NumUsefulEvts << " matched with HiSCORE events" << endl;

    DrawSepEvt(evts, 8074, 1);
    FindDeadTime(evts, 1);
    return 0;
}