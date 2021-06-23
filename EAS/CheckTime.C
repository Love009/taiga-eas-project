#include "common.h"
#include "ReadTreeOfEvts.h"
#include <iostream>
#include <dirent.h>
#include <sys/types.h>
using namespace std;

static const int    N_STN_TOT =  19;
static const int    N_STN_MIN =   2;
static const int    DT_WINDOW = 50000; // ns
const char *        DATA_DIRECTORY         ="/home/liuba/workspace/TAIGA/taiga-eas/EAS/data_Grande+rex/271119.03.rsg/";
const char*         HiSCORE_DATA_ADDRESS   = "HiSCORE_data/2019-2020/out_2611.dat";
//===============================================================================================//
//===============================================================================================//
//===============================================================================================//
int HiscTime(const char* hiscDat = HiSCORE_DATA_ADDRESS)
{
    vector<Event> evts, evtsNew;
    string line, lineLast;
    ifstream hiscore(hiscDat);
    if (! hiscore.is_open()) return -1;
    int count = 0;
    int T1[2], T2[2];
    const char u = ':';
    while(getline(hiscore, line))
    {
        count++;
        int k = 0;
        while(char(line[k]) != u)  k++;
        if(count == 1)
        {
            T1[0] = stoi(line.substr(k - 2, 9));
            T1[1] = stoi(line.substr(k + 1, 2));
        }
        lineLast = line;
    }
    int k = 0;
    while(char(lineLast[k]) != u)  k++;
    T2[0] = stoi(lineLast.substr(k - 2, 9));
    T2[1] = stoi(lineLast.substr(k + 1, 2));
    hiscore.close();
    cout << "HiSCORE time interval: " << T1[0] << ":" << T1[1] << " - " << T2[0] << ":" << T2[1] << endl;
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Read binary file (data from station № stID) and add events to vector evts
int ReadBinary(string fileName, int stID, vector<Event> &evts, int fileNum)
{
#include "grandeFormat.h"
    Station st;
    Event evt;
    int Trow[4];
    int stT = 0;
    ifstream myfile(fileName, std::ios::binary);
    if(!myfile){
        if(fileNum == 1) cout << "File doesn't exist:" << fileName << endl;
        return -1;
    }
    else{
        while(! myfile.eof()){
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
            if(h < 24 && m < 60 && s < 60 && mks < 1000 && mls < 1000)
            {
                Trow[0] = h;
                Trow[1] = m;
                Trow[2] = s;
                Trow[3] = ns + 1000 * mks + 1000000 * mls;
                stT = ns + 1000 * mks + 1000000 * mls;
            }
            else
            {
                cout << "Time format error   " << Trow[0] << ":" << Trow[1] << ":" << Trow[2] << ":" << mls << ":" << mks << ":" << ns << endl;
                break;
            }
            //Data reading
            for(int u = 0; u < 8; u++){
                myfile.read(reinterpret_cast<char*>(&VMEadress), sizeof(UInt_t));
                myfile.read(reinterpret_cast<char*>(&ndata), sizeof(UShort_t));
                myfile.read(reinterpret_cast<char*>(&channel), sizeof(UShort_t));
                data = (UShort_t *) realloc (data, ndata * sizeof(UShort_t ));
                myfile.read(reinterpret_cast<char*>(data), ndata * sizeof(UShort_t));
            }
            myfile.read(reinterpret_cast<char*>(&end), 4);//End of package    FFFFFFFF 
            if(end != 0xffffffff)    cout<<"Data error in file:  " << fileName << endl;
            else if(evtID == evtIDlast){}
            //Writing to events vector
            else{
                st.SetID(stID);
                st.SetT(stT+Delay);
                st.SetTrow(Trow);
                st.SetDelay(Delay);
                evt.Reset();
                evt.AddStation(st);
                Trow[3] += Delay;
                evt.SetEvtTime(Trow);
                evtIDlast = evtID;
                evts.push_back(evt);
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
    if (dir == NULL) return;
    while ((entry = readdir(dir)) != NULL) fileNames.push_back(entry->d_name);
    closedir(dir);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Build a vector of events from Grande
void BuildEvts(vector<Event> &evts, const char * dir)
{
    int stID = 1;
    string fileName, temp;
    temp = string(dir) + "St31";
    vector<string> fileNames;
    list_dir(temp.c_str(), fileNames);
    string suffix = "/" + fileNames[(int)fileNames.size() - 1].substr(0, 6);
    while(stID <= N_STN_TOT){
        int folderNum = stID + 30;
        int fileNum = 1, end_ = 0;
        do{
            fileName =  string(dir) + "St" + to_string(folderNum) + suffix + to_string(folderNum);
            if(fileNum < 10) fileName += ".00"+ to_string(fileNum);
            else if(fileNum < 100) fileName += ".0"+ to_string(fileNum);
            else    fileName += "."+ to_string(fileNum);
            end_ = ReadBinary(fileName, stID, evts, fileNum);
            fileNum++;
        }while(end_ == 0);
        stID += N_STN_TOT - 1;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Write tree of events coinciding with HiSCORE to root-file
int TimeInterval(const char * dir = DATA_DIRECTORY, const char * hiscAdress = HiSCORE_DATA_ADDRESS)
{
    gROOT->ProcessLine("#include <vector>");
    vector<Event> evts;
    BuildEvts(evts, dir);
    cout << "\t\t\t   hh:mm\nGrande time interval : " << evts[0].MuStn(0).Trow()[0] << ":" << evts[0].MuStn(0).Trow()[1] << " - " << evts[(int)evts.size() - 1].MuStn(0).Trow()[0] << ":" <<  evts[(int)evts.size() - 1].MuStn(0).Trow()[1] << endl;
    HiscTime(hiscAdress);
    return 0;
}
