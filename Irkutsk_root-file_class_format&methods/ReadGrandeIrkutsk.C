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

vector<double> MuonX[2];
vector<double> MuonY[2];
vector<double> MuonZ[2];
vector<double> MuonTerr;
int trigNum = 0, NumWrite = 0, LastNumWrite = 0;//numbers of packages
//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
//Initualize stations coordinates and time errors (vectors MuonX, Y, Z, Terr)
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
        }
        else{
            MuonX[1].push_back(x);
            MuonY[1].push_back(y);
            MuonZ[1].push_back(z);
        }
        i++;
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CheckWordErr(UShort_t *data, UShort_t n)//res = 0 - OK, res = 1 - error
{
    int res = 0;
    for(int i = 0; i < n; i++)
    {
        UShort_t var;
        for(int l = 1; l <= 9; l++)
        {
            var = 0xfef0|l;
            if(data[i] == var) {cout << "word error in " << i << " word "; res = 1; break;}
        }
        var = 0xfef0|0xa;
        if(!(data[i] & var)) {cout << "word error in " << i << " word "; res = 1; break;}
        var = 0xfef0|0xb;
        if(!(data[i] & var)) {cout << "word error in " << i << " word "; res = 1; break;}
        var = 0xfef0|0xc;
        if(!(data[i] & var)) {cout << "word error in " << i << " word "; res = 1; break;}
        var = 0xfef0|0xd;
        if(!(data[i] & var)) {cout << "word error in " << i << " word "; res = 1; break;}
        var = 0xfef0|0xe;
        if(!(data[i] & var)) {cout << "word error in " << i << " word "; res = 1; break;}
        var = 0xfef0|0xf;
        if(!(data[i] & var)) {cout << "word error in " << i << " word "; res = 1; break;}
    }
    return res;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Add station st to existing evt in vector evts or add new evt with station st (all these by comparing trigger times)
int WriteCloseEvt(vector<Event> &evts, Station st, int stID)
{
    long int i = 0;
    long int dT;
    do
    {
        long int dT0 = (evts[i].MuStn(0).Traw()[0] - st.Traw()[0]) * 3600 +(evts[i].MuStn(0).Traw()[1] - st.Traw()[1]) * 60 + (evts[i].MuStn(0).Traw()[2] - st.Traw()[2]);
        dT = dT0 * 1000000000 + evts[i].MuStn(0).Traw()[3] + evts[i].MuStn(0).Delay() - (st.Traw()[3] + st.Delay());
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
        evts.insert(evts.begin() + i, ev);            
        NumWrite++;
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Read binary file (data from station № stID) and add events to vector evts
int ReadBinary(string fileName, int stID, vector<Event> &evts, string * date)
{
    #include "grandeFormat.h"
    Station st;
    Event evt;
    vector<int> chan[ADCchanNum];
    int Traw[4];
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
            Delay = ((clNumber>>16)&0xffff ) * 10/ 2.; // cable length in ns
           // cout<<Delay<<endl;
            if(h < 24 && m < 60 && s < 60 && mks < 1000 && mls < 1000)
            {
                Traw[0] = h;
                Traw[1] = m;
                Traw[2] = s;
                Traw[3] = ns + 1000 * mks + 1000000 * mls;
            }
            else
            {
                cout << "Time format error   " << Traw[0] << ":" << Traw[1] << ":" << Traw[2] << ":" << mls << ":" << mks << ":" << ns << endl;
                break;
            }
    //Data reading
            for(int j = 0; j < ADCchanNum; j++) {chan[j].clear(); ADCwordError[j] = 1;}
            for(int u = 0; u < ADCchanNum; u++){
                myfile.read(reinterpret_cast<char*>(&VMEadress), sizeof(UInt_t));
                myfile.read(reinterpret_cast<char*>(&ndata), sizeof(UShort_t));
                myfile.read(reinterpret_cast<char*>(&channel), sizeof(UShort_t));
                data = (UShort_t *) realloc (data, ndata * sizeof(UShort_t ));
                myfile.read(reinterpret_cast<char*>(data), ndata * sizeof(UShort_t));
    //Saving of signal from trigger channels
                int wordErr = CheckWordErr(data, ndata);
                if(evtID != evtIDlast)
                { 
                    for(int i = 0; i < (int) ndata; i++)  chan[u].push_back((int)data[i]);//What about errors in words???
                    if(wordErr) { ADCwordError[u] = 0;  cout << "ADC channel " << u << ": word error " ;for(int k=0; k < ADCchanNum; k++)cout << ADCwordError[k];cout<< "\n" << endl;}
                }
            }
            myfile.read(reinterpret_cast<char*>(&end), 4);//End of package    FFFFFFFF 
            if(end != 0xffffffff)    cout<< "Data error in file:  " << fileName << endl;
            else if(evtID == evtIDlast){}
    //Writing to events vector
            else{
                trigNum++;
                double x[2], y[2], z[2];
                x[0] = MuonX[0][stID - 1];
                x[1] = MuonX[1][stID - 1];
                st.SetX(x);
                y[0] = MuonY[0][stID - 1];
                y[1] = MuonY[1][stID - 1];
                st.SetY(y);
                z[0] = MuonZ[0][stID - 1];
                z[1] = MuonZ[1][stID - 1];
                st.SetZ(z);

                st.SetID(stID);
                st.SetTraw(Traw);
                st.SetDelay(Delay);
                st.ResetADC();
                for(int ii = 0; ii < ADCchanNum; ii++)  st.SetADC(ii, chan[ii]);
                st.SetADCwordErr(ADCwordError);
                evt.Reset();
                evt.AddStation(st);
                evtIDlast = evtID;
                if(stID == 1)  {
                    evt.SetDate(date);
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
    //gInterpreter->GenerateDictionary("Event","TunkaIrkutsk.h");
    int stID, folderNum = 31;
    string fileName, temp;
    Initialize();
    trigNum = 0;
    NumWrite = 0;
    
    temp = string(dir) + "St31";
    vector<string> fileNames;
    list_dir(temp.c_str(), fileNames);
    string suffix = "/" + fileNames[(int)fileNames.size() - 1].substr(0, 6);
    string date[3];//yy:mm:dd
    date[2] = fileNames[(int)fileNames.size() - 1].substr(0, 2);//day
    date[1] = fileNames[(int)fileNames.size() - 1].substr(2, 2);//month
    date[0] = fileNames[(int)fileNames.size() - 1].substr(4, 1);//year
    if(date[0] == '0') date[0] = '2' + date[0];
    else date[0] = '1' + date[0];
    do{
        stID = folderNum - 30;
        int fileNum = 1, end_ = 0;
        do{
            fileName =  string(dir) + "St" + to_string(folderNum) + suffix + to_string(folderNum);
            if(fileNum < 10) fileName += ".00"+ to_string(fileNum);
            else if(fileNum < 100) fileName += ".0"+ to_string(fileNum);
            else    fileName += "."+ to_string(fileNum);
            end_ = ReadBinary(fileName, stID, evts, date);
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
    }
    cout << "evts SIZE_after N_STN_MIN selection: " << evts.size() << endl;
    int newNumWrite = 0;
    for(int i = 0; i < (int)evts.size(); i++) newNumWrite += evts[i].NStations();
    cout << "N_packs_read = " << trigNum << "\t N_packs_written = " << NumWrite << "\tN_packs_rest (N_trig_stations >= N_STN_MIN)= " << newNumWrite<< endl;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
    cout << evts[0].MuStn(0).Traw()[0] << " -  " << evts[(int)evts.size() - 1].MuStn(0).Traw()[0] << " hours" << endl;
    for(int i = 0; i < (int)evts.size(); i++)
    {
        event = & evts[i];
        Muon.Fill();
    }
    f.Write();
    f.Close();
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
