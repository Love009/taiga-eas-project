#include "common.h"
#include <map>
int plata1 = 33, plata2 = 62;//ns
int N_CLASTERS = 3;

int JoinTXTfiles()
{
    ofstream out("XYZDelaysHiSCORE.txt");
    ifstream inXYZ("coordinatesHiSCORE.txt");
    ifstream inDelay1("delaysHiSCORE.txt");//individual
    ifstream inDelay2("delay_between_HiSC.dat");//betw clasters and data collection center
   // ifstream inDelay3("gro10_status.txt");//electronic type delay

    string title;
    double numLen;
    double id1, id2, id3, x, y, z, delay1, delay2[3], delay3;
    int DRS;
    getline(inXYZ, title);
    out << "    " << title << "	" << "Delay, ns" << "\n";
    getline(inDelay1, title);
    getline(inDelay2, title);
    //getline(inDelay3, title);
//Read optic length to data collection center
    for(int u = 0; u < N_CLASTERS; u++)
    {
        double rest;
        inDelay2 >> rest >> delay2[u] >> rest;
        cout << delay2[u] << endl;
    }
    inDelay2.close();

//Read coords and two other types of delay and write to the joined file
    while (1) {
    //coords
        inXYZ >> id1 >> x >> y >> z;
    //optic length to claster center
        do{
            inDelay1 >> id2 >> delay1 >> numLen;
            cout << id2 << "\t" << delay1 << endl;
        }while(id2 != id1);
 /*   //electronic board delay
        do{
            inDelay3 >> id3 >> DRS;
        }while(id3 != id1);
        if(DRS < 100.) delay3 = plata1;
        else delay3 = plata2;*/
    //write to file
        if (!inXYZ.good()) break;

        if(numLen != 0) {
            int clasterNumber = int(id1 / 100.);
            out << id1 << "\t";
            out << /*fixed << setprecision(2) <<*/ x << "\t" << y << "\t" << z << "\t"
                << delay1 + delay2[clasterNumber] /*+ delay3*/ << "\n";
        }
    }
    inDelay1.close();
    //inDelay3.close();
    out.close();
    cout << "FINISH" << endl;
    return 0;
}