#include <map>
#include "common.h"

vector<float> HiscX;
vector<float> HiscY;
vector<float> HiscZ;
vector<float> delayTot;
map <int, double[4]> HiSCcoordsAndDelays;

void InitualizeStCoordsAndDelays()
{
    double x, y, z, t, delay;
    vector<double> Delay;
    string  title, id;
    int i = 1;

    ifstream in1("XYZDelaysHiSCORE.txt");
    getline(in1, title);
    while (1) {
        in1 >> id >> x >> y >> z;
        if (!in1.good()) break;
        double XYZDelay[4];
        XYZDelay[0] = x;
        XYZDelay[1] = y;
        XYZDelay[2] = z;
        XYZDelay[3] = delay;
        i++;
    }

    HiSCcoordsAndDelays.insert (pair<int, double[4]> (id, XYZDelay));
    for (auto it = HiSCcoordsAndDelays.begin(); it != HiSCcoordsAndDelays.end(); ++it)
    {
        cout << (*it).first << " : " << (*it).second[0] << endl;
    }
    return 0;
}

int CorrectDelays()
{

    return 0;
}