#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;
#include "TMath.h"

int Calculate(const char *FileName)
{
    int n3 = 0, n4 = 0, nHiSC3 = 0, nHiSC4 = 0, flagHiSCend = 0;
    double tGr1[4], tGr2[4], tHiSC1[4], tHiSC2[4];
    vector<Event> evts;
    ReadTreeOfEvts(FileName, evts);
    cout << evts.size() << endl;
    for(int l = 0; l < 4; l++) tGr1[l] = evts[0].MuStn(0).Trow()[l];
    for(int l = 0; l < 4; l++) tGr2[l] = evts[(int)evts.size() - 1].MuStn(0).Trow()[l];
    for(int i = 0; i < (int)evts.size(); i++)
    {
        if(evts[i].NStations() >= 2)
        {
            n3++;
            if(evts[i].MatchHiSC() == 1)
            {
                if(nHiSC3 == 0)
                {
                    for(int j = 0; j < 4; j++) tHiSC1[j] = evts[i].EvtTimeHiSC()[j];
                }
                else flagHiSCend = i;
                nHiSC3++;
            }
        }
        if(evts[i].NStations() >= 4)
        {
            n4++;
            if(evts[i].MatchHiSC() == 1) nHiSC4++;
        }
    }
    for(int j = 0; j < 4; j++) tHiSC2[j] = evts[flagHiSCend].EvtTimeHiSC()[j];
    evts.clear();

    double dtGr = (tGr2[0] - tGr1[0]) * 60 * 60 + (tGr2[1] - tGr1[1]) * 60 + (tGr2[2] - tGr1[2]);
    dtGr += dtGr + (tGr2[3] - tGr1[3]) / 1000000000.;
    dtGr = dtGr / 3600.;

    double dtHiSC = (tHiSC2[0] - tHiSC1[0]) * 60 * 60 + (tHiSC2[1] - tHiSC1[1]) * 60 + (tHiSC2[2] - tHiSC1[2]);
    dtHiSC += dtHiSC + (tHiSC2[3] - tHiSC1[3]) / 1000000000.;
    dtHiSC = dtHiSC / 3600.;

    double GrRate3 = n3 / dtGr;
    double GrRate4 = n4 / dtGr;
    double GrHiscRate = nHiSC4 / dtHiSC;
    cout << "Gr rate3 = " << GrRate3 << " evts/hour\tGr rate4 = " << GrRate4 << " evts/hour\n";
    cout << "gr-hisc rate = " << GrHiscRate << " evts/hour\n";

    return 0;
}