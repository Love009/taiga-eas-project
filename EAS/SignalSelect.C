//CorrectedTree(const char * FileName = "analysis_results/Gr+Hisc271119.root", const char * FileToWrite = "analysis_results/test.root", int threshold = 5, int lowLimit = -10) is the main function
//which writes a tree from the file FileName with correced time to the file FileToWrite
#include "ReadTreeOfEvts.h"
using namespace std;

const int TIME_WINDOW = 5; // *5 ns
const int LOW_SIG_CHECK_INTERVAL = 30; // counts of 200 total ADS counts
static const int    N_STN_MIN =   2;
vector<double> MeanCH0;
vector<double> MeanCH2;
vector<double> SigmaCH0;
vector<double> SigmaCH2;
const char *INITIALIZATION_FILE_NAME = "analysis_results/pedestals/pedestal261119.txt";
const char *INPUT_FILE_NAME          = "analysis_results/Gr+Hisc291219.root";
const char *OUTPUT_FILE_NAME         = "analysis_results/test291219.root";
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
            //TCanvas *c1 = new TCanvas("c1", "Signal from channel");
            //TH1I *sig = new TH1I("sig","Signal from channel 0, too low", (int)ADCsig.size(), 0, (int)ADCsig.size());
            //for(int k = 0; k < (int)ADCsig.size(); k++) { sig->SetBinContent((k + 1), ADCsig[k]-MeanCH0[stID-1]);}
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
int CorrectedTree(const char * InitFileName = INITIALIZATION_FILE_NAME, const char * FileName = INPUT_FILE_NAME, const char * FileToWrite = OUTPUT_FILE_NAME, int threshold = 5, int lowLimit = -10)
{
    gErrorIgnoreLevel=kError;
    vector<Event> evts;
    ReadTreeOfEvts(FileName, evts);
    Initialize(InitFileName);
    cout << SigmaCH0[0] << endl;
    int nSigTot = 0, nSig = 0, sigTime0/*ns*/, sigTime2/*ns*/, counter = 0, counter_init = 0, evtsBefore = evts.size();
    for(int i = 0; i < (int) evts.size(); i++)
    {
       for(int j = 0; j < evts[i].NStations(); j++)
        {
           int stID = evts[i].MuStn(j).ID();
           counter_init++;
         
           sigTime0 = FindSignal(evts[i].MuStn(j).ADC0(), stID, 0, SigmaCH0[stID - 1], threshold, lowLimit);
            sigTime2 = FindSignal(evts[i].MuStn(j).ADC2(), stID, 2, SigmaCH2[stID - 1], threshold, lowLimit);
             if(sigTime0 > 0 && sigTime2 > 0)
              {
                  int sigTime = sigTime0;
                  if(sigTime2 < sigTime) sigTime = sigTime2;
                  double newT = evts[i].MuStn(j).T() + (double)sigTime;
                  evts[i].MuStn(j).SetT(newT);
                  counter++;

                  //TCanvas *c3 = new TCanvas("c3", "Signal from channel");
                  //TH1I *sig3 = new TH1I("sig3","Signal from channel 0, signal IS found", (int)evts[i].MuStn(0).ADC0().size(), 0, (int)evts[i].MuStn(0).ADC0().size());
                  //cout << "Station " << stID << " Signal from channel 0, signal is found"<< endl;
                  //for(int k = 0; k < (int) evts[i].MuStn(j).ADC0().size(); k++) { sig3->SetBinContent((k + 1), evts[i].MuStn(j).ADC0()[k]-MeanCH0[stID-1]);}
                  //sig3 -> Draw("");
              }
              else
              {
                  if(sigTime0 < 0 && sigTime0 != -2*5)
                  {
                     // TCanvas *c0 = new TCanvas("c0", "Signal from channel");
                      TH1I *sig0 = new TH1I("sig0","Signal from channel 0, signal IS NOT found", (int)evts[i].MuStn(0).ADC0().size(), 0, (int)evts[i].MuStn(0).ADC0().size());
                      for(int k = 0; k < (int) evts[i].MuStn(j).ADC0().size(); k++) { sig0->SetBinContent((k + 1), evts[i].MuStn(j).ADC0()[k]-MeanCH0[stID-1]);} //cout << evts[i].MuStn(j).ADC2()[k]-MeanCH2[stID-1]<<endl;}
                      //sig0 -> Draw("");
                  }
                  if(sigTime2 < 0 && sigTime2 != -2*5)
                  {
                      //TCanvas *c2 = new TCanvas("c2", "Signal from channel");
                      TH1I *sig2 = new TH1I("sig2","Signal from channel 2, signal IS NOT found", (int)evts[i].MuStn(0).ADC2().size(), 0, (int)evts[i].MuStn(0).ADC2().size());
                      for(int k = 0; k < (int) evts[i].MuStn(j).ADC2().size(); k++) { sig2->SetBinContent((k + 1), evts[i].MuStn(j).ADC2()[k]-MeanCH2[stID-1]); } //cout << evts[i].MuStn(j).ADC2()[k]-MeanCH2[stID-1]<<endl;}
                      //sig2 -> Draw("");
                  }

                  if(evts[i].NStations() <= N_STN_MIN)
                  {
                      cout << i << " event is ERASED\tsigTime0 = " << sigTime0 << "\tsigTime2 = " << sigTime2 << "\tStation_number = " << stID << "\tNstns = " << evts[i].NStations() << endl;
                      evts.erase(evts.begin() + i);
                      i--;
                      break;
                  }
                  else
                  {
                      evts[i].DeleteStation(j);//delete j-th station from the i-th event
                      cout << "DELETE STATION " << j + 1 << " from event " << i << endl;
                      j--;
                  }

              }
        }
    }
    cout << "Selected " << counter << " of " << counter_init << " signals" << endl;
    cout << evtsBefore << " events before\t" << evts.size() << " events after"<< endl;
//Write new root-file   
    TFile f(FileToWrite,"recreate");
    TTree Muon("Muon","");
    gROOT->ProcessLine("#include <vector>");
    Event *event = new Event();
    Muon.Branch("events", &event);
    Muon.SetAutoSave(1000000);
    for(int i = 0; i < (int)evts.size(); i++)
    {
        event = & evts[i];
        Muon.Fill(); 
    }
    f.Write();
    f.Close();
    return 0;    
}
