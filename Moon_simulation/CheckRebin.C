#include "common.h"
#include "ReadTreeOfEvts.h"
using namespace std;
#include "TMath.h"
int n_bins = 800.;

//auto c52 = new TCanvas("c52", "");

double CalcMean(TH2D *hist, int n)
{
    double mean = 0;
    int counter = 0;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            double tmp0 = hist -> GetBinContent(i, j);
            if(tmp0 != 0){ mean += tmp0; counter++;}
        }
    }
    mean = mean / counter;
    return mean;
}
//-------------------------------------------------------------------------------------
double CalcRMS(TH2D *hist, double mean, int n)
{
    double rms = 0;
    int counter = 0;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            double tmp0 = hist -> GetBinContent(i, j);
            if(tmp0 != 0) {rms += (tmp0 - mean) * (tmp0 - mean); counter++;}
        }
    }
    rms = sqrt(rms / counter);
    return rms;
}
//-------------------------------------------------------------------------------------
int ShiftHistContent(TH2D *hist, TH2D *hist2, double mean)
{
    int n = hist->GetXaxis()->GetNbins();
    cout << "n in shift = " << n << endl;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            int bin = hist -> GetBin(i, j);
            double tmp = (hist -> GetBinContent(bin)) - mean;
            hist2 -> SetBinContent(bin, tmp);
        }
    }
    return 0;
}
//-------------------------------------------------------------------------------------
int CheckBinning(TH2D *hist, double mean, int n)
{
    int flag = 1;
    double more = 0, less = 0;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            int bin = hist -> GetBin(i, j);
            double tmp = (hist -> GetBinContent(bin)) - mean;
            if(tmp < 0) less++;
            else more++;
        }
    }
    cout << "l/m = " << less / more << endl;
    if(less / more >= 0.8 && less / more <= 1.2) flag = 0;
    cout << "less = " << less << "  more = " << more << endl;
    return flag;
}
//-------------------------------------------------------------------------------------
int ScaleToSigma(TH2D *theta_phi_gen, TH2D *theta_phi_rec)
{
    int flag1 = 0, flag2 = 0;
    int rebin_num = 0;
    double mean_rec, mean_gen, rms_rec, rms_gen;
    int n;
    do {
        if(flag1 == 1 || flag2 == 1) {
            theta_phi_gen->Rebin2D(2,2);
            theta_phi_rec->Rebin2D(2,2);
            rebin_num++;
            flag1 = 0;
            flag2 = 0;
        }
        n = theta_phi_gen->GetXaxis()->GetNbins();
        cout <<"rebin = "<< rebin_num << "  n=" <<n <<endl;
        mean_gen = CalcMean(theta_phi_gen, n);
        mean_rec = CalcMean(theta_phi_rec, n);
        rms_gen = CalcRMS(theta_phi_gen, mean_gen, n);
        rms_rec = CalcRMS(theta_phi_rec, mean_rec, n);
        if(rebin_num <= 1) flag1 = CheckBinning(theta_phi_gen, mean_gen, n);
        if(rebin_num <= 1) flag2 = CheckBinning(theta_phi_rec, mean_rec, n);
        cout << "mean_gen = " << mean_gen << "  rms_gen = " << rms_gen << endl;
        cout << "mean_rec = " << mean_rec << "  rms_rec = " << rms_rec << "\n" <<endl;
    }while(flag1 == 1 || flag2 == 1);

    n = theta_phi_gen->GetXaxis()->GetNbins();
    TH2D *theta_phi_rec2 = new TH2D("theta_phi_rec2", "theta_phi_rec2", n, -6., 6., n, -6., 6.);
    TH2D *theta_phi_gen2 = new TH2D("theta_phi_gen2", "theta_phi_gen2", n, -6., 6., n, -6., 6.);
    auto c53 = new TCanvas("c53", "");
    c53->Divide(2,1);
    ShiftHistContent(theta_phi_gen, theta_phi_gen2, mean_gen);
    ShiftHistContent(theta_phi_rec, theta_phi_rec2, mean_rec);
    Double_t scale_rec = 1. / rms_rec;
    Double_t scale_gen = 1. / rms_gen;
    theta_phi_rec2->Scale(scale_rec);
    theta_phi_gen2->Scale(scale_gen);
    c53->cd(1);
    theta_phi_rec2 -> Draw("CONT4Z");
    c53->cd(2);
    theta_phi_gen2 -> Draw("CONT4Z");
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main1()
{
    TH2D *theta_phi_rec = new TH2D("theta_phi_rec", "theta_phi_rec", n_bins, -6., 6., n_bins, -6., 6.);
    TH2D *theta_phi_gen = new TH2D("theta_phi_gen", "theta_phi_gen", n_bins, -6., 6., n_bins, -6., 6.);
    TRandom2 *r1=new TRandom2();
    TRandom2 *r2=new TRandom2();
    for(int i = 0; i < 600000; i++)
    {
        double valx = (gRandom->Rndm()) * n_bins;
        double valy = (gRandom->Rndm()) * n_bins;
        double valz = 0.;
        if(gRandom->Rndm() > 1/2. && (valx - n_bins / 2.)*(valx-n_bins / 2.) + (valy-n_bins / 2.)*(valy-n_bins / 2.) >= 0.3 * n_bins) valz = 40. + 10. * r1->Gaus(0.,1.);
        theta_phi_rec -> SetBinContent(valx, valy, valz);
    }
    for(int i = 0; i < 600000; i++)
    {
        double valx = (gRandom->Rndm()) * n_bins;
        double valy = (gRandom->Rndm()) * n_bins;
        double valz = 0.;
        if(gRandom->Rndm() > 1/2. && (valx - n_bins / 2.)*(valx-n_bins / 2.) + (valy-n_bins / 2.)*(valy-n_bins / 2.) >= 0.3 * n_bins) valz = 40. + 2 * r1->Gaus(0.,1.);
        theta_phi_gen -> SetBinContent(valx, valy, valz);
    }
  /*  auto c52 = new TCanvas("c52", "");
    c52->Divide(2,1);
    c52->cd(1);
    theta_phi_rec -> GetXaxis() -> SetTitle("(PhiRec - PhiMoon)*Sin(ThetaMoon)");
    theta_phi_rec -> GetYaxis() -> SetTitle("ThetaRec - ThetaMoon");
    theta_phi_rec -> Draw("CONT4Z");
    c52->cd(2);
    theta_phi_gen -> GetXaxis() -> SetTitle("(PhiGen - PhiMoon)*Sin(ThetaMoon)");
    theta_phi_gen -> GetYaxis() -> SetTitle("ThetaGen - ThetaMoon");
    theta_phi_gen -> Draw("CONT4Z");*/
    ScaleToSigma(theta_phi_gen, theta_phi_rec);
    return 0;
}