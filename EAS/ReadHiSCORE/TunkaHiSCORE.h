//gInterpreter->GenerateDictionary("Event","Tunka.h");
//gInterpreter->GenerateDictionary("Station","Tunka.h");
#include <stdlib.h>
#include<iostream>
#include<vector>
using namespace std;

class StationHiSCORE
{
private:
    int   sID;
    double sX;
    double sY;
    double sZ;
    double sT;//stT = ns + 1000 * mks + 1000000 * mls + Delay// sTrow[3] only
    double sDelay[3];//ns, optic length to claster center, optic length to center of data collection, electronic delay
    int sTrow[4];//hh:mm:ss,ns

public:
    Station()
    {
        sID = -1;
        sX = - 1000000.;
        sY = - 1000000.;
        sZ = - 1000000.;
        sT = - 1000000.;
        sDelay[0] = 0.;
        sDelay[1] = 0.;
        sDelay[2] = 0.;
        sTrow[0] = -1;
        sTrow[1] = -1;
        sTrow[2] = -1;
        sTrow[3] = -1;
    }

    int     ID()        {  return sID;     }
    double  X()         {  return  sX;     }
    double  Y()         {  return  sY;     }
    double  Z()         {  return  sZ;     }
    double  T()         {  return  sT;     }
    double *Delay()     {  return sDelay;  }
    int    *Trow()	    {  return sTrow;   }

    void SetID(int i)       { sID = i;         }
    void SetX(double t)     {  sX = t;         }
    void SetY(double t)     {  sY = t;         }
    void SetZ(double t)     {  sZ = t;         }
    void SetT(double t)     {  sT = t;         }
    void SetDelay(double t[3]) {for(int i = 0; i < 3; i++)  sDelay[i] = t[i]; }
    void SetTrow(int t[4]){for(int i = 0; i < 4; i++) sTrow[i] = t[i];}
};

class EventHiSCORE
{
private:

    double the_gen;
    double phi_gen;
    double the_rec;
    double phi_rec;
    double chi2;
    int type;
    int matchHiSC;
    int evtTime[4];
    int evtTimeHiSC[4];
    vector<Station> Stn;
    double sEAScenterRec[2];//(x,y) coords
    double sEAScenterGen[2];

public:
    Event()
    {
        matchHiSC      = 0;  //0-no info from HiSC, 1 - there is HiSC info
        evtTime[0] = -1;
        evtTime[1] = -1;
        evtTime[2] = -1;
        evtTime[3] = -1;
        evtTimeHiSC[0] = -1;
        evtTimeHiSC[1] = -1;
        evtTimeHiSC[2] = -1;
        evtTimeHiSC[3] = -1;
    }
    // Member Functions()

    int      Type()                  {  return type;       }
    int      NStations()             {  return Stn.size(); }
    int      *EvtTime()              {  return evtTime;    }
    int      *EvtTimeHiSC()          {  return evtTimeHiSC;}
    int      MatchHiSC()             {  return matchHiSC;    }
    Station& MuStn(int i)            {  return Stn[i];     }

    double   Chi2()                  {  return chi2;       }
    double   ThetaRec()              {  return the_rec;    }
    double   ThetaGen()              {  return the_gen;    }
    double   PhiRec()                {  return phi_rec;    }
    double   PhiGen()                {  return phi_gen;    }
    double  *EAScentersRec()       {  return sEAScenterRec;}
    double  *EAScentersGen()       {  return sEAScenterGen;}


    void     AddStation(Station &s, int i = -1)  { if(i == -1)Stn.push_back(s);
        else Stn.insert(Stn.begin() + i, s);}
    void     DeleteStation(int i){Stn.erase(Stn.begin() + i);}
    void     SetType(int i)          { type    = i;        }
    void     SetChi2(double t)       { chi2    = t;        }
    void     SetThetaRec(double t)   { the_rec = t;        }
    void     SetThetaGen(double t)   { the_gen = t;        }
    void     SetPhiRec(double t)     { phi_rec = t;        }
    void     SetPhiGen(double t)     { phi_gen = t;        }
    void     SetEvtTime(int t[4])    {for(int i = 0; i < 4; i++) evtTime[i] = t[i];    }
    void     SetEvtTimeHiSC(int t[4]){for(int i = 0; i < 4; i++) evtTimeHiSC[i] = t[i];}
    void     SetEAScenterRec(double t[]) {sEAScenterRec[0] = t[0]; sEAScenterRec[1] = t[1];}
    void     SetEAScenterGen(double t[]) {sEAScenterGen[0] = t[0]; sEAScenterGen[1] = t[1];}
    void     SetMatchHiSC(int i)     { matchHiSC     = i;        }

    void     Reset()     { Stn.clear();        }
};

,