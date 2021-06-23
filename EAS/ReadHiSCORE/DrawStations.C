TH2D *coords_gr = new TH2D("coords_gr", "coords_gr", 2 * 900, -900, 900, 2 * 900, -900, 900);
TH2D *coords_hisc = new TH2D("coords_hisc", "coords_hisc", 2 * 900, -900, 900, 2 * 900, -900, 900);
TH2I *st_number = new TH2I("st_number", "number of station", 2 * 900, -900, 900, 2 * 900, -900, 900);

void InitializeGrande(int par = 1)
{
    double x, y, z, sigma;
    string  title, rest;
    int i = 1;
    ifstream in("coordinatesWithTimeErr.txt");
    getline(in, title);
    int count = 0;
    while (1) {
        in >> rest >> x >> y >> z >> sigma;// sigma = t_error
        if (!in.good()) break;
        if(i % 2 == par % 2){//2%2 for getting underground counters coordinstes, 1%2 for onground
            coords_gr -> Fill(y, x);
            //cout << x << "  " << y << "  " << count + 1 << endl;
            st_number -> SetBinContent(y + 900., x + 900., count + 1);
           // cout << count << ":" << x <<"  " << y << endl;
            count++;
        }
        i++;
    }
}

void InitializeHiSCORE()
{
    double x, y, z;
    string  title, rest;
    ifstream in("coordinatesHiSCORE.txt");
    getline(in, title);
    int count = 0;
    while (1) {
        in >> rest >> x >> y >> z;
       // cout << rest << endl;
        st_number -> SetBinContent(-y + 900., x + 900, stoi(rest));
        if (!in.good()) break;
        coords_hisc -> Fill(-y, x);
        //cout << count << ":" << x <<"  " << y << endl;
        count++;
    }
}

//==============================================================================
int DrawSts()
{
    InitializeGrande();
    InitializeHiSCORE();
    coords_gr -> SetMarkerSize(1);
    coords_hisc -> SetMarkerSize(1);
    auto c2 = new TCanvas("c2", "");
    coords_gr -> SetMarkerStyle(21);
    coords_gr -> SetMarkerColorAlpha(kPink, 0.35);
    coords_gr -> Draw();
    coords_hisc -> SetMarkerStyle(21);
    coords_hisc -> SetMarkerColorAlpha(kGreen, 0.35);
    coords_hisc -> Draw("SAME");
    st_number -> SetMarkerSize(1.2);
    //st_number -> SetMarkerColor(kRed);
    st_number -> Draw("SAME TEXT");
    return 0;
}