int WriteCoordWithTimeErr(){//to the .txt file
    double sigma = 15;//ns 
    double x, y, z;
    string  title, rest;
    ofstream out("txt-files/coordinatesWithTimeErr.txt");
    ifstream in("coordinatesGRANDE.txt"); 
  
    getline(in, title);
    out << "      	" << title << "	" << "sigma, ns" << "\n";
    while (1) {
        in >> rest >> x >> y >> z;
        if (!in.good()) break;
        out << rest << "	";
        out << fixed << setprecision(2) << x << "	" << y << "	" << z << "	" << sigma << "\n";
    }
    in.close();
    out.close(); 
    return 0;
}
