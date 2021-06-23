//gInterpreter->GenerateDictionary("Event","TunkaIrkutsk.h");
#include <stdlib.h>
#include<iostream>
#include<vector>
using namespace std;
static const int ADCchanNum = 8;

class Station
{
private:
    int   sID;
    double sX;
    double sY;
    double sZ;
    int sDelay;   //ns
    int sTraw[4];//hh:mm:ss,ns
    vector<int> sADC[ADCchanNum];//static array of vectors
    int sADCwordError[ADCchanNum];//sADCwordError[0]<->ADC0, sADCwordError[2]<->ADC2... sADCwordError[i] = 0 => Word error in ADCi trace

public:
    Station()
    {
        sID = -1;
        sX = - 1000000.;
        sY = - 1000000.;
        sZ = - 1000000.;
        sDelay = 0;//ns
        sTraw[0] = -1;
        sTraw[1] = -1;
        sTraw[2] = -1;
        sTraw[3] = -1; //Treal = Traw + Delay
        for(int i = 0; i< ADCchanNum; i++)sADCwordError[i] = 1;

    }

    int ID()                {  return sID;     }
    double X()              {  return  sX;     }
    double Y()              {  return  sY;     }
    double Z()              {  return  sZ;     }
    int Delay()             {  return sDelay;  }
    int * Traw()	        {  return sTraw;   }
    vector<int> ADC(int i)  {  return sADC[i]; }
    vector<int> *ADC()      {  return sADC;    }
    int *ADCwordErr()       {  return sADCwordError;}

    void SetID(int i)       { sID = i;       }
    void SetX(double t)     {  sX = t;       }
    void SetY(double t)     {  sY = t;       }
    void SetZ(double t)     {  sZ = t;       }
    void SetDelay(int t)    { sDelay = t;    }
    void SetTraw(int t[4])  {for(int i = 0; i < 4; i++) sTraw[i] = t[i];}
    void SetADC(int chan, vector<int> t) {for(int i = 0; i < (int) t.size(); i++) sADC[chan].push_back(t[i]);}
    void SetADCwordErr(int t[ADCchanNum]){for(int i = 0; i < ADCchanNum; i++)sADCwordError[i] = t[i];}
    void ResetADC()              {for(int i = 0; i < ADCchanNum; i++) sADC[i].clear(); }

};


class Event
{
private:

    string date[3];      //yy:mm:dd
    int matchHiSC;       //0-no info from HiSC, 1 - there is HiSC info
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
        matchHiSC          = 0  ;
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
    int      NStations()             {  return Stn.size(); }
    string   *Date()                 {  return date;       }
    double   *XYHiSC()               {  return xy_HiSC;    }
    int      *EvtTimeHiSC()          {  return evtTimeHiSC;}
    Station& MuStn(int i)            {  return Stn[i];     }
    double   ThetaHiSC()             {  return the_HiSC;   }
    double   PhiHiSC()               {  return phi_HiSC;   }
    double   RaHiSC()                {  return ra_HiSC;    }
    double   DecHiSC()               {  return dec_HiSC;   }
    int      MatchHiSC()             {  return matchHiSC;  }

    void     AddStation(Station &s, int i = -1)  { if(i == -1)Stn.push_back(s);
                                                   else Stn.insert(Stn.begin() + i, s);}
    void     DeleteStation(int i)    {Stn.erase(Stn.begin() + i);}
    void     SetMatchHiSC(int i)     { matchHiSC    = i;         }
    void     SetThetaHiSC(double t)  { the_HiSC = t;             }
    void     SetPhiHiSC(double t)    { phi_HiSC = t;             }
    void     SetRaHiSC(double t)     { ra_HiSC = t;              }
    void     SetDecHiSC(double t)    { dec_HiSC = t;             }
    void     SetEvtTimeHiSC(int t[4]){ for(int i = 0; i < 4; i++) evtTimeHiSC[i] = t[i];}
    void     SetDate(string t[3])    { for(int i = 0; i < 3; i++) date[i] = t[i];       }
    void     SetXYHiSC(double t[2])  { for(int i = 0; i < 1; i++) xy_HiSC[i] = t[i];    }
    void     Reset()                 { Stn.clear();              }
};