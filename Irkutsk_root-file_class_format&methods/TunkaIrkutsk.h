//gInterpreter->GenerateDictionary("Event","TunkaIrkutsk.h");
#include <stdlib.h>
#include<iostream>
#include<vector>
#include "common.h"
using namespace std;
static const int ADCchanNum    = 8 ;
static const int ADCpointsNum  = 20; //For pedestal parameter calculation
static const int N_STN_TOT     = 19;
const int LOW_SIG_CHECK_INTERVAL = 30; // counts of 200 total ADS counts

class Station
{
private:
    int   sID;
    double sX[2];//ground and underground
    double sY[2];
    double sZ[2];
    int sDelay;   //ns
    int sTraw[4];//hh:mm:ss,ns
    vector<int> sADC[ADCchanNum];//static array of vectors
    int sADCwordError[ADCchanNum];//sADCwordError[0]<->ADC0, sADCwordError[2]<->ADC2... sADCwordError[i] = 0 => Word error in ADCi trace

public:
    Station()
    {
        sID = -1;
        sX[0] = - 1000000.;
        sY[0] = - 1000000.;
        sZ[0] = - 1000000.;
        sX[1] = - 1000000.;
        sY[1] = - 1000000.;
        sZ[1] = - 1000000.;
        sDelay = 0;//ns
        sTraw[0] = -1;
        sTraw[1] = -1;
        sTraw[2] = -1;
        sTraw[3] = -1; //Treal = Traw + Delay
        for(int i = 0; i< ADCchanNum; i++)sADCwordError[i] = 1;

    }

    int ID()                {  return sID;     }
    double *X()             {  return  sX;     }
    double *Y()             {  return  sY;     }
    double *Z()             {  return  sZ;     }
    int Delay()             {  return sDelay;  }
    int *Traw()	            {  return sTraw;   }
    vector<int> ADC(int i)  {  return sADC[i]; }
    vector<int> *ADC()      {  return sADC;    }
    int *ADCwordErr()       {  return sADCwordError;}

    void SetID(int i)       {  sID = i;       }
    void SetX(double *t)    {  sX[0] = t[0]; sX[1] = t[1];}
    void SetY(double *t)    {  sY[0] = t[0]; sY[1] = t[1];}
    void SetZ(double *t)    {  sZ[0] = t[0]; sZ[1] = t[1];}
    void SetDelay(int t)    { sDelay = t;    }
    void SetTraw(int t[4])  {for(int i = 0; i < 4; i++) sTraw[i] = t[i];}
    void SetADC(int chan, vector<int> t) {for(int i = 0; i < (int) t.size(); i++) sADC[chan].push_back(t[i]);}
    void SetADCwordErr(int t[ADCchanNum]){for(int i = 0; i < ADCchanNum; i++)sADCwordError[i] = t[i];}
    void ResetADC()         {for(int i = 0; i < ADCchanNum; i++) sADC[i].clear(); }

};


class Event
{
private:

    string date[3];      //yy:mm:dd
    int matchHiSC;       //0-no info from HiSC, 1 - there is HiSC info
    int matchTunka133;
    int matchIACT;
    double the_HiSC;
    double phi_HiSC;
    double ra_HiSC;     //sky coords: right ascension,
    double dec_HiSC;    //declination
    int evtTimeHiSC[4];
    double xy_HiSC[2];  //EAS center coordinates
    vector<Station> Stn;

public:
    Event()
    {
        matchHiSC      = 0      ;
        matchTunka133  = 0      ;
        matchIACT      = 0      ;
        the_HiSC       = - 1000.;
        phi_HiSC       = - 1000.;
        date[0]        = "-1"   ;
        date[1]        = "-1"   ;
        date[2]        = "-1"   ;
        evtTimeHiSC[0] = -1     ;
        evtTimeHiSC[1] = -1     ;
        evtTimeHiSC[2] = -1     ;
        evtTimeHiSC[3] = -1     ;
        xy_HiSC[0]     = - 1000.;
        xy_HiSC[1]     = - 1000.;
    }

    // Member Functions()
    int      NStations()             {  return Stn.size();   }
    string   *Date()                 {  return date;         }
    double   *XYHiSC()               {  return xy_HiSC;      }
    int      *EvtTimeHiSC()          {  return evtTimeHiSC;  }
    Station& MuStn(int i)            {  return Stn[i];       }
    double   ThetaHiSC()             {  return the_HiSC;     }
    double   PhiHiSC()               {  return phi_HiSC;     }
    double   RaHiSC()                {  return ra_HiSC;      }
    double   DecHiSC()               {  return dec_HiSC;     }
    int      MatchHiSC()             {  return matchHiSC;    }
    int      MatchTunka133()         {  return matchTunka133;}
    int      MatchIACT()             {  return matchIACT;    }

    void     AddStation(Station &s, int i = -1)  { if(i == -1)Stn.push_back(s);
                                                   else Stn.insert(Stn.begin() + i, s);}
    void     DeleteStation(int i)    {Stn.erase(Stn.begin() + i);}
    void     SetMatchHiSC(int i)     { matchHiSC     = i;        }
    void     SetMatchTunka133(int i) { matchTunka133 = i;        }
    void     SetMatchIACT(int i)     { matchIACT     = i;        }
    void     SetThetaHiSC(double t)  { the_HiSC = t;             }
    void     SetPhiHiSC(double t)    { phi_HiSC = t;             }
    void     SetRaHiSC(double t)     { ra_HiSC = t;              }
    void     SetDecHiSC(double t)    { dec_HiSC = t;             }
    void     SetEvtTimeHiSC(int t[4]){ for(int i = 0; i < 4; i++) evtTimeHiSC[i] = t[i];}
    void     SetDate(string t[3])    { for(int i = 0; i < 3; i++) date[i] = t[i];       }
    void     SetXYHiSC(double t[2])  { for(int i = 0; i < 2; i++) xy_HiSC[i] = t[i];    }
    void     Reset()                 { Stn.clear();              }
};

//Extended class
//need a method for pedestal parameter calculation
class EventAnalysis: public Event
{
private:
    double the_Gr;
    double phi_Gr;
    double chi2;
    long int evtTime[4];
    double xy_Gr[2];  //EAS center coordinates obtained from Grande data
    //vector<vector<double>> &evtPedestalMean(ADCchanNum , vector<double> (N_STN_TOT));
    double evtPedestalMean[N_STN_TOT][ADCchanNum];//for each channel of each station
    double evtPedestalRMS[N_STN_TOT][ADCchanNum];
    int sigStart[N_STN_TOT][ADCchanNum];//ns (-1000 => no sig)
    int sigMax[N_STN_TOT][ADCchanNum];
public:
    EventAnalysis()
    {
        the_Gr = -1000.;
        phi_Gr = -1000.;
        chi2 = -10;
        evtTime[0] = -1;
        evtTime[1] = -1;
        evtTime[2] = -1;
        evtTime[3] = -1;
        xy_Gr[0]     = - 1000.;
        xy_Gr[1]     = - 1000.;
    }
    EventAnalysis(Event &evt)
    {
        the_Gr = -1000.;
        phi_Gr = -1000.;
        chi2 = -10;
        evtTime[0] = -1;
        evtTime[1] = -1;
        evtTime[2] = -1;
        evtTime[3] = -1;
        xy_Gr[0]   = - 1000.;
        xy_Gr[1]   = - 1000.;

        for(int i = 0; i < evt.NStations(); i++)
        {
            this->AddStation(evt.MuStn(i));
        }
        this->SetDate(evt.Date());
        this->SetMatchHiSC(evt.MatchHiSC());
        this->SetEvtTimeHiSC(evt.EvtTimeHiSC());
        this->SetXYHiSC(evt.XYHiSC());
        this->SetPhiHiSC(evt.PhiHiSC());
        this->SetThetaHiSC(evt.ThetaHiSC());
        this->SetRaHiSC(evt.RaHiSC());
        this->SetDecHiSC(evt.DecHiSC());
    }

    // Member Functions()
    long int *EvtTime()                 {  return evtTime; }
    double   Chi2()                     {  return chi2;    }
    double   ThetaGr()                  {  return the_Gr;  }
    double   PhiGr()                    {  return phi_Gr;  }
    double   *XYGr()                    {  return xy_Gr;   }
    double   *EvtPedestalMean(int stNum){  return evtPedestalMean[stNum];}//0<=stNum<=18
    double   *EvtPedestalRMS(int stNum) {  return evtPedestalRMS[stNum]; }
    int      *SigStart(int stNum)       {  return sigStart[stNum];       }
    int      *SigMax(int stNum)            {  return sigMax[stNum];         }

    void     SetEvtPedestalMean(double t[N_STN_TOT][ADCchanNum]){for(int i=0; i < N_STN_TOT; i++)for(int j=0; j < ADCchanNum; j++)evtPedestalMean[i][j]=t[i][j];}
    void     SetEvtPedestalRMS(double t[N_STN_TOT][ADCchanNum]){for(int i=0; i < N_STN_TOT; i++)for(int j=0; j < ADCchanNum; j++)evtPedestalRMS[i][j]=t[i][j];}
    void     SetEvtPedestalMean(double t[N_STN_TOT], int ADCnum){for(int i=0; i < N_STN_TOT; i++)evtPedestalMean[i][ADCnum]=t[i];}
    void     SetEvtPedestalRMS(double t[N_STN_TOT], int ADCnum){for(int i=0; i < N_STN_TOT; i++)evtPedestalRMS[i][ADCnum]=t[i];}
    void     SetChi2(double t)        { chi2   = t;      }
    void     SetThetaGr(double t)     { the_Gr = t;      }
    void     SetPhiGr(double t)       { phi_Gr = t;      }
    void     SetEvtTime(long int t[4]){ for(int i = 0; i < 4; i++) evtTime[i] = t[i]; }
    void     SetXYGr(double t[2])     { for(int i = 0; i < 2; i++) xy_Gr[i] = t[i];   }
    void     CalculatePedestalPars()  {

        double mean[N_STN_TOT][ADCchanNum], RMS[N_STN_TOT][ADCchanNum];
        int norm[N_STN_TOT][ADCchanNum];
        for(int k = 0; k < N_STN_TOT; k++) for(int l = 0; l < ADCchanNum; l++){
            mean[k][l] = 0.;
            RMS[k][l] = 0.;
            norm[k][l] = ADCpointsNum;
        }
        int Nsts = this->NStations();

        for(int j = 0; j < Nsts; j++){
            int stID = this->MuStn(j).ID() - 1;
            Station& st = this->MuStn(j);
            for(int l = 0; l < ADCchanNum; l++){
                for(int k = 0; k < ADCpointsNum; k++){
                        if(st.ADC(l)[k] < 2048 + 50 && st.ADC(l)[k] > 2048 - 50) mean[stID][l] += st.ADC(l)[k];
                        else norm[stID][l]--;
                }
                evtPedestalMean[stID][l] = mean[stID][l] / norm[stID][l];
            }   
        }

        for(int k = 0; k < N_STN_TOT; k++) for(int l = 0; l < ADCchanNum; l++) norm[k][l] = ADCpointsNum;

        for(int j = 0; j < Nsts; j++){
            int stID = this->MuStn(j).ID() - 1;
            Station& st = this->MuStn(j);
            for(int l = 0; l < ADCchanNum; l++){
                for(int k = 0; k < ADCpointsNum; k++){
                        if(st.ADC(l)[k] < 2048 + 50 && st.ADC(l)[k] > 2048 - 50) RMS[stID][l] += (st.ADC(l)[k] - evtPedestalMean[stID][l]) * (st.ADC(l)[k] - evtPedestalMean[stID][l]);
                        else norm[stID][l]--;
                }
                evtPedestalRMS[stID][l] = sqrt(RMS[stID][l] / norm[stID][l]);
            }   
        }
    }
    void FindSignal(int threshold = 5 /*sigmas*/, int lowLimit = -10 /*sigmas*/, int window = 5/*counts*/){
            int nSts = this->NStations();
            for(int stNum = 0; stNum < nSts; stNum++)
            {
                Station st = this->MuStn(stNum);
                int stID = st.ID() - 1;
                for(int chanNum = 0; chanNum < ADCchanNum; chanNum++)
                {
                    int sigTime = -200, flag = 0;
                    double mean = evtPedestalMean[stID][chanNum];
                    double sigma = evtPedestalRMS[stID][chanNum];
                    for(int k = 0; k < (int) st.ADC(chanNum).size(); k++)
                    {
                        int delta = st.ADC(chanNum)[k] -  mean;
                        if(delta <= lowLimit * sigma && k <= LOW_SIG_CHECK_INTERVAL) 
                        {
                           // cout << "Signal ampl is lower than " << lowLimit * sigma << "\tStation_number = " << stID << "\t channel " << chanNum << "\n\n";

                            //TH1I *sig = new TH1I("sig","Signal from channel 0, too low", (int)st.ADC(chanNum).size(), 0, (int)st.ADC(chanNum).size());
                            //for(int k = 0; k < (int)st.ADC(chanNum).size(); k++) { sig->SetBinContent((k + 1), delta);}
                            // delete sig;
                            sigTime = -200;
                            break;
                        }
                        else if(flag >= window)
                        {
                            TArrayD x(2*window), y(window);
                            TMatrixD A(window, 2), A_t(2, window),Y(window, 1), Pars(2, 1);
                            for(int i = 0; i < window; i++)
                            {
                                y[i] = st.ADC(chanNum)[k - 1 - window + i] - mean;
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
                    sigStart[stID][chanNum] = sigTime * 5;
                }
            }
        }
        void FindSigMax()
        {
            int oscSize = (int)this->MuStn(0).ADC(0).size();
            TH1I *osc = new TH1I("osc", "oscillogramm", oscSize, 0, oscSize);
            for(int i = 0; i < this-> NStations(); i++)
            {
                int stationID = this->MuStn(i).ID() - 1;
                for(int l = 0; l < ADCchanNum; l++)
                {
                    osc -> Reset();
                    for (int j = 0; j < oscSize; j++) osc->SetBinContent(j + 1, this->MuStn(i).ADC(l)[j]);
                    int maxValue = osc->GetMaximum();
                    if(sigStart[stationID][l] != -1000. && maxValue < 4096-10) sigMax[stationID][l] = maxValue - this -> EvtPedestalMean(stationID)[l];
                    else sigMax[stationID][l] = -1000.;
                }

            }
            delete osc;
        }

};
