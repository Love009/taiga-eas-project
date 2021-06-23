#include "common.h"
#include "ReadTreeOfEvts.h"
#include <iostream>
#include <dirent.h>
#include <sys/types.h>
using namespace std;

static const int    N_STN_TOT =  19;
static const int    N_STN_MIN =   2;
static const int    DT_WINDOW = 10000; // ns
const char *        DATA_DIRECTORY               ="/home/liuba/workspace/TAIGA/taiga-eas/EAS/data_Grande+rex/271119.03.rsg/";
const char *        OUTPUT_FILE_NAME            = "analysis_results/test.root";

vector<float> MuonX;
vector<float> MuonY;
vector<float> MuonZ;
vector<float> MuonTerr;
int trigNum = 0, NumWrite = 0, LastNumWrite = 0;//numbers of packages

TH2F *deltaT = new TH2F("deltaT", "deltaT, ns", 1000., -10., 3000., 1000., -10., 3000.);

//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
//Initualize stations coordinates and time errors (vectors MuonX, Y, Z, Terr)
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//For culc of random coincidences
int less1mks, more1mks;
void FindTimeMIN_MAX(Event &evt, float *res)
{
    float *dtime = new float[evt.NStations() - 1];
    float t[N_STN_TOT];
    for(int i = 0; i < evt.NStations(); i++)
    {
        cout << evt.MuStn(i).Trow()[0] << ":" << evt.MuStn(i).Trow()[1] << ":" << evt.MuStn(i).Trow()[2] << ":"<<evt.MuStn(i).Trow()[3] << endl;
        t[i] = (float)evt.MuStn(i).Trow()[3];
    }
    sort(t, t + evt.NStations());
    for(int i = evt.NStations() - 1; i > 0; i--) dtime[i - 1] = t[i] - t[i - 1];
    sort(dtime, dtime + evt.NStations() - 1);

    res[0] = dtime[0];
    res[1] = dtime[evt.NStations() - 2];

    if(res[1] >= 1500.) more1mks++;
    else less1mks++;
    for(int l = 0; l < evt.NStations() - 1; l++)
    {
        cout << dtime[l] << "\t";
    }
    cout << "\n";
    delete [] dtime;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Add station st to existing evt in vector evts or add new evt with station st (all these by comparing trigger times)
int WriteCloseEvt(vector<Event> &evts, Station st, int stID)
{
    long int i = 0;
    long int dT;
    do
    {
        long int dT0 = (evts[i].EvtTime()[0] - st.Trow()[0]) * 3600 +(evts[i].EvtTime()[1] - st.Trow()[1]) * 60 + (evts[i].EvtTime()[2] - st.Trow()[2]);
        dT = dT0 * 1000000000 + evts[i].EvtTime()[3] - (st.Trow()[3] + st.Delay());
        int flagg = 0;
        for(int u = 0; u < evts[i].NStations(); u++)  if(evts[i].MuStn(u).ID() == st.ID()) flagg = 1;
        if(abs(dT) <= DT_WINDOW /*ns*/&& flagg == 0) 
        {
            evts[i].AddStation(st);
            NumWrite++;
            break;
        }
        else if(dT > 0)
        { 
            Event ev;
            ev.AddStation(st);
            int *T = st.Trow();
            T[3] += st.Delay();
            ev.SetEvtTime(T);
            evts.insert(evts.begin() + i, ev);
            NumWrite++;
            break;
        }
        else i++;
    }while(i < (long int) evts.size());
    if(LastNumWrite == NumWrite)//This station trigger time is more than all time values recorded in the vector
    {
        Event ev;
        ev.AddStation(st);
        int *T = st.Trow();
        T[3] += st.Delay();
        ev.SetEvtTime(T);
        evts.insert(evts.begin() + i, ev);            
        NumWrite++;
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Read binary file (data from station № stID) and add events to vector evts
int ReadBinary(string fileName, int stID, vector<Event> &evts)
{
    #include "grandeFormat.h"
    Station st;
    Event evt;
    vector<int> chan0, chan2, chan1, chan3;
    int Trow[4]; 
    int stT = 0;
    ifstream myfile(fileName, std::ios::binary);  
    if(! myfile){
        cout << "End of the folder:" << fileName << endl;
        return -1;
    }
    else{ 
        while(! myfile.eof()){        
            nPacksTot++;
            myfile.read(reinterpret_cast<char*>(&linkID), sizeof(UShort_t));
            myfile.read(reinterpret_cast<char*>(&packSize), sizeof(UShort_t));
            myfile.read(reinterpret_cast<char*>(&evtID), sizeof(UInt_t));
            myfile.read(reinterpret_cast<char*>(&errorsTot), sizeof(UInt_t));
            myfile.read(reinterpret_cast<char*>(&time), sizeof(ULong_t ));
            myfile.read(reinterpret_cast<char*>(&clNumber), sizeof(UInt_t));
    //Time reading    hh:mm:ss,mls.mks.ns               
            for(int i = 0; i < 4; i++)  data_time[i]=((time >> (48 - 16*(3 - i))) & 0xffff);
            ns = ( data_time[0] & 0x7f ) * 10; // 7 bits for ns
            mks = ( data_time[0] & 0xff80 ) >> 7;
            mks |= (data_time[1] & 1) << 9;// (9 + 1) bits for mks
            mls = ( data_time[1] & 0x7fe ) >> 1;// 10 bits for mls
            s   = ( data_time[1] & 0xf800 ) >> 11;
            s |= (data_time[2] & 1) << 5;// (5 + 1) bits for seconds
            m  = ( data_time[2] & 0x7e ) >> 1;// 6 bits for minutes
            h = ( data_time[2] & 0xf80 ) >> 7; // 5 bits for hours
            Delay = ((clNumber>>16)&0xffff ) * 10 / 2.; // cable length in ns 
            //cout <<"ID = " << clNumber << "\tRead delay = " << Delay << endl;
            if(h < 24 && m < 60 && s < 60 && mks < 1000 && mls < 1000)
            {
                Trow[0] = h;
                Trow[1] = m;
                Trow[2] = s;
                Trow[3] = ns + 1000 * mks + 1000000 * mls;
               // cout << h << ":" << m << ":" << s << ":" << ns << endl;
                stT = ns + 1000 * mks + 1000000 * mls;
            }
            else
            {
                cout << "Time format error   " << Trow[0] << ":" << Trow[1] << ":" << Trow[2] << ":" << mls << ":" << mks << ":" << ns << endl;
                break;
            }
    //Data reading
            chan0.clear();//chan1*10
            chan2.clear();//chan3*10
            chan1.clear();
            chan3.clear();
            for(int u = 0; u < 8; u++){  
                myfile.read(reinterpret_cast<char*>(&VMEadress), sizeof(UInt_t));
                myfile.read(reinterpret_cast<char*>(&ndata), sizeof(UShort_t));
                myfile.read(reinterpret_cast<char*>(&channel), sizeof(UShort_t));
                data = (UShort_t *) realloc (data, ndata * sizeof(UShort_t ));
                myfile.read(reinterpret_cast<char*>(data), ndata * sizeof(UShort_t));
    //Saving of signal from trigger channels   
                if((u <= 3)&& evtID != evtIDlast)
                { 
                    for(int i = 0; i < (int) ndata; i++) 
                    {
                        if(u == 0) chan0.push_back((int)data[i]);
                        else if(u == 1) chan1.push_back((int)data[i]);
                        else if(u == 2) chan2.push_back((int)data[i]);
                        else chan3.push_back((int)data[i]);
                    }
                }
            }
            myfile.read(reinterpret_cast<char*>(&end), 4);//End of package    FFFFFFFF 
            if(end != 0xffffffff)    cout<< "Data error in file:  " << fileName << endl;
            else if(evtID == evtIDlast){}
    //Writing to events vector
            else{
                trigNum++;
                st.SetX(MuonX[stID - 1]);
                st.SetY(MuonY[stID - 1]);
                st.SetZ(MuonZ[stID - 1]);
                st.SetID(stID);
                st.SetT(stT+Delay);
                st.SetTerr( MuonTerr[stID - 1]);
                st.SetTrow(Trow);
                st.SetDelay(Delay);
                
                st.ResetADC();
                st.SetADC0(chan0);
                st.SetADC2(chan2);
                st.SetADC1(chan1);
                st.SetADC3(chan3);
                
                evt.Reset();
                evt.AddStation(st);
    //Drawing the signal of one of the packs
               /*  if(NumWrite == 1)
                {
                    TCanvas *c1 = new TCanvas("c1", "Signal from channeles");
                    c1 -> Divide(2,1);
                    TH1I *sig0 = new TH1I("sig0","Signal from channel 0", ndata, 0, ndata);
                    TH1I *sig2 = new TH1I("sig2","Signal from channel 2", ndata, 0, ndata);
                    for(int i = 0; i < (int) evt.MuStn(0).ADC0().size(); i++) sig0 -> SetBinContent((i + 1), evt.MuStn(0).ADC0()[i]);
                    for(int i = 0; i < (int) evt.MuStn(0).ADC2().size(); i++) sig2 -> SetBinContent((i + 1), evt.MuStn(0).ADC2()[i]);
                    c1 -> cd(1);
                    sig0 -> Draw();
                    c1 -> cd(2);
                    sig2 -> Draw();
                }*/
                Trow[3] += Delay;
                evt.SetEvtTime(Trow);
                evtIDlast = evtID;
                if(stID == 1)  {
                    evts.push_back(evt);
                    NumWrite++;
                }
                else  { 
                    LastNumWrite = NumWrite; 
                    WriteCloseEvt(evts, st, stID);
                }
            }
        }
        myfile.close();
        return 0;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Reads file names from the directory path and fills fileNames array
void list_dir(const char *path, vector<string> &fileNames) {
    struct dirent *entry;
    DIR *dir = opendir(path);

    if (dir == NULL) {
        return;
    }
    while ((entry = readdir(dir)) != NULL) fileNames.push_back(entry->d_name);
    closedir(dir);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Build a vector of >= N_STN_MIN station events from Grande
void BuildEvts(vector<Event> &evts, const char * dir)
{
    //gInterpreter->GenerateDictionary("Event","Tunka.h");    
    int stID, folderNum = 31;
    string fileName, temp;
    Initialize(1);// par = 1 <- onground stations coordinates
    trigNum = 0;
    NumWrite = 0;
    
    temp = string(dir) + "St31";
    vector<string> fileNames;
    list_dir(temp.c_str(), fileNames);
    string suffix = "/" + fileNames[(int)fileNames.size() - 1].substr(0, 6);
    do{
        stID = folderNum - 30;
        int fileNum = 1, end_ = 0;
        do{
            fileName =  string(dir) + "St" + to_string(folderNum) + suffix + to_string(folderNum);
            if(fileNum < 10) fileName += ".00"+ to_string(fileNum);
            else if(fileNum < 100) fileName += ".0"+ to_string(fileNum);
            else    fileName += "."+ to_string(fileNum);
            end_ = ReadBinary(fileName, stID, evts);
            fileNum++;
        }while(end_ == 0);
        folderNum++;
    }while(stID < N_STN_TOT);

    cout << "evts SIZE_before: " << evts.size() << endl;
    for(int i = 0; i < (int) evts.size(); i++) 
    {
        if(evts[i].NStations() < N_STN_MIN) {
            evts.erase(evts.begin() + i);
            i--;
        }
       /* else
        {
            for(int kk = 0; kk < evts[i].NStations(); kk++)
            {
                cout << "\t" << i << ":   " << evts[i].MuStn(kk).ID() << endl;
                cout << evts[i].MuStn(kk).Trow()[0] << ":" << evts[i].MuStn(kk).Trow()[1]  << ":" << evts[i].MuStn(kk).Trow()[2]  << ":" << evts[i].MuStn(kk).Trow()[3]  << endl;
            }
        }*/
    }
    cout << "evts SIZE_after N_STN_MIN selection: " << evts.size() << endl;
    cout << "evts SIZE_after N_STN_MIN selection: " << evts.size() << endl;
    int newNumWrite = 0;
    for(int i = 0; i < (int)evts.size(); i++) newNumWrite += evts[i].NStations();
    cout << "N_packs_read = " << trigNum << "\t N_packs_written = " << NumWrite << "\tN_packs_rest (N_trig_stations >= N_STN_MIN)= " << newNumWrite<< endl;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Write tree of N_STN_MIN events coinciding with HiSCORE to root-file
int ReadGrande(const char * dir = DATA_DIRECTORY, const char *  fileToWrite = OUTPUT_FILE_NAME)
{
    gROOT->ProcessLine("#include <vector>");
    vector<Event> evts;
    BuildEvts(evts, dir);
    TFile f(fileToWrite,"recreate");
    TTree Muon("Muon","");
    gROOT->ProcessLine("#include <vector>");
    Event *event = new Event();
    Muon.Branch("events", &event);
    Muon.SetAutoSave(1000000);
    cout << evts[0].MuStn(0).Trow()[0] << " -  " << evts[(int)evts.size() - 1].MuStn(0).Trow()[0] << " hours" << endl;

    less1mks = 0;
    more1mks = 0;
    for(int i = 0; i < (int)evts.size(); i++)
    {
        if(evts[i].NStations() == 3) {
            float timeMIN_MAX[2];
            FindTimeMIN_MAX(evts[i], timeMIN_MAX);
            //cout << timeMIN_MAX[0]<< "\t" << timeMIN_MAX[1] << endl;
            deltaT->Fill(timeMIN_MAX[0], timeMIN_MAX[1]);
        }
        event = & evts[i];
        Muon.Fill();
    }
    f.Write();
    f.Close();

    auto c100 = new TCanvas("c100","");
    deltaT -> Draw();
    cout << "more1mks part = " << more1mks / ((less1mks + more1mks) * 1.) * 100. << " percent" << endl;
    return 0; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void JoinTrees(string FileName1 = "3StEvts.root", string FileName2 = "3StEvts1.root", string FileToWrite = "27.11.19.3StGrandeEvts+HiSCORE.root")
{
    vector<Event> evts;
    ReadTreeOfEvts(FileName1.c_str(), evts);
    ReadTreeOfEvts(FileName2.c_str(), evts);
    
    TFile f(FileToWrite.c_str(),"recreate");
    TTree Muon("Muon","");
    gROOT->ProcessLine("#include <vector>");
    Event *event = new Event();
    Muon.Branch("events", &event);
    Muon.SetAutoSave(1000000);
    for(int i = 0; i <(int) evts.size(); i++)
    {
        event = & evts[i];
        Muon.Fill(); 
    }     
    f.Write();
    f.Close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
    ReadGrande("data_Grande+rex/010120.01.rsg/", "analysis_results/Gr010120.01.root");
    ReadGrande("data_Grande+rex/010120.02.rsg/", "analysis_results/Gr010120.02.root");

    JoinTrees("analysis_results/Gr010120.01.root", "analysis_results/Gr010120.02.root", "analysis_results/Gr010120.01-02.root");
    return 0;
}
