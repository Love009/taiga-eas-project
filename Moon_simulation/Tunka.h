//gInterpreter->GenerateDictionary("Event","Tunka.h"); 
//gInterpreter->GenerateDictionary("Station","Tunka.h");  
#include <stdlib.h>
#include<iostream>
#include<vector>
using namespace std;

class Station
{

private:
  int   sID;
  double sX;
  double sY;
  double sZ;
  double sT;//stT = ns + 1000 * mks + 1000000 * mls + Delay// sTrow[3] only
  double sTerr;
  double sDelay;//ns
  int sTrow[4];//hh:mm:ss,ns
  vector<int> sADC0;
  vector<int> sADC2;
  vector<int> sADC1;
  vector<int> sADC3;
  //ADCdat *sE = new ADCdat();
  
public: 
  Station()
  {
        sID = -1;
        sX = - 1000000.;
        sY = - 1000000.;
        sZ = - 1000000.;
        sT = - 1000000.;
        sTerr = - 1000000.;
        sDelay = 0.;
        sTrow[0] = -1;
        sTrow[1] = -1;
        sTrow[2] = -1;
        sTrow[3] = -1;
  }
  
  int ID()                {  return sID;     }
  double X()              {  return  sX;     }
  double Y()              {  return  sY;     }
  double Z()              {  return  sZ;     }
  double T()              {  return  sT;     }
  double Terr()           {  return  sTerr;  }
  double Delay()          {  return sDelay;  }
  int * Trow()	          {  return sTrow;   }
  vector<int> ADC0(){return sADC0;}
  vector<int> ADC2(){return sADC2;}
  vector<int> ADC1(){return sADC1;}
  vector<int> ADC3(){return sADC3;}
  //ADCdat * E() {return sE;}
  
  void SetID(int i)       { sID = i;         }
  void SetX(double t)     {  sX = t;         }
  void SetY(double t)     {  sY = t;         }
  void SetZ(double t)     {  sZ = t;         }
  void SetT(double t)     {  sT = t;         }
  void SetTerr(double t)  {  sTerr = t;      }
  void SetDelay(double t){ sDelay = t;    }
  void SetTrow(int t[4]){for(int i = 0; i < 4; i++) sTrow[i] = t[i];}
  void SetADC0(vector<int> t) {for(int i = 0; i < (int) t.size(); i++) sADC0.push_back(t[i]);}
  void SetADC2(vector<int> t) {for(int i = 0; i < (int) t.size(); i++) sADC2.push_back(t[i]);}
  void SetADC1(vector<int> t) {for(int i = 0; i < (int) t.size(); i++) sADC1.push_back(t[i]);}
  void SetADC3(vector<int> t) {for(int i = 0; i < (int) t.size(); i++) sADC3.push_back(t[i]);}
  //void SetE(ADCdat *t){for(int i = 0; i < (int) t->Chan0.size(); i++) sE->Chan0.push_back(t->Chan0[i]);
  //for(int i = 0; i < (int) t->Chan2.size(); i++) sE->Chan2.push_back(t->Chan2[i]);}
  void ResetADC()             {sADC0.clear(); sADC2.clear(); sADC1.clear(); sADC3.clear();}
  
};

class Event
{ 
private:

  double the_gen;
  double phi_gen;  
  double the_rec;
  double phi_rec;
  double chi2;  
  int type;  
  int evtTime[4];
  int evtTimeHiSC[4];
  vector<Station> Stn;
  double sEAScenterRec[2];//(x,y) coords
  double sEAScenterGen[2];

public: 
  Event()
  {
      the_gen = - 10.;
      phi_gen = -10.;
      the_rec = -10.;
      phi_rec = -10.;
      chi2 = -10;
      type = 0;
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

  void     Reset()     { Stn.clear();        }
}; 

