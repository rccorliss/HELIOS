//
// Template functions to make nice looking plot
//
// Axel Drees 1/1/2015
//            9/17/2021 converted to class structure
//

#include "MyPlot.h"


MyPlot::MyPlot(){ }               // default constructor
MyPlot::~MyPlot(void){}

TCanvas *MyPlot::Canvas(const Char_t *name, Int_t px2=720, Int_t py2 =920, Int_t px1=10,Int_t py1=10, Int_t logy=0, Int_t logx=0 ){
//
// prepare canvas with default options
    TCanvas *cCanvas = new TCanvas(name,"",px1,py1,px2,py2); // set plot size and aspect ratio
    cCanvas->SetTickx(1);                                  // ticks on top
    cCanvas->SetTicky(1);                                  // ticks on right
    cCanvas->SetBorderSize(0);                             // border size
    cCanvas->SetFillColor(0);
    cCanvas->SetMargin(LeftMargin,RightMargin,BottomMargin,TopMargin);                 // set margins
    if(logy) cCanvas->SetLogy();                           // 0 lin, 1 log
    if(logx) cCanvas->SetLogx();  
    return cCanvas;
}

TH1D *MyPlot::Frame( const Char_t *name, const Char_t *xaxis, const Char_t *yaxis, Double_t xmin=0., Double_t xmax=10., Double_t ymin=0., Double_t ymax=10., Int_t center=0){
//
 // setup empty histogram as frame
    TH1D *hFrame1D = new TH1D(name,"", 100, xmin, xmax);      // define x axis of plot
    hFrame1D->SetMinimum(ymin);                               // minimum y
    hFrame1D->SetMaximum(ymax);                               // maximum y
    hFrame1D->SetStats(false);
//
// setup y axis of plot 
    hFrame1D->GetYaxis()->SetTitleFont(Font);
    hFrame1D->GetYaxis()->SetTitleSize(TitleSize);                       
    hFrame1D->GetYaxis()->SetTitleOffset(yTitleOffset);
    if (center) hFrame1D->GetYaxis()->CenterTitle();
    hFrame1D->GetYaxis()->SetLabelFont(Font);
    hFrame1D->GetYaxis()->SetLabelSize(LabelSize);
    hFrame1D->GetYaxis()->SetTitle(yaxis);
//    
// setup x axis
    hFrame1D->GetXaxis()->SetTitleFont(Font);
    hFrame1D->GetXaxis()->SetTitleSize(TitleSize);
    hFrame1D->GetXaxis()->SetTitleOffset(xTitleOffset);
    if (center) hFrame1D->GetXaxis()->CenterTitle();
    hFrame1D->GetXaxis()->SetLabelFont(Font);
    hFrame1D->GetXaxis()->SetLabelSize(LabelSize); 
    hFrame1D->GetXaxis()->SetTitle(xaxis);
    return hFrame1D;
}
//

TLegend *MyPlot::Legend(const Char_t *name, Double_t x1=0.6, Double_t y1=0.79,Double_t x2=0.8,Double_t y2=0.94,Double_t size=1.){
//
// prepare legend
//
  TLegend *lLegend = new TLegend(x1,y1,x2,y2,name);   // define and locate (relative coordinate 0-1)
  lLegend->SetBorderSize(0);                          // no border
  lLegend->SetFillColor(0);                           // no color
  lLegend->SetTextSize(LegendSize*size);                     // set ledgend size as needed 
  lLegend->SetTextFont(Font);                           // font size should be consistent with PlotFrame
  lLegend->SetTextColor(LegendColor);
  return lLegend;
}

void MyPlot::StyleMe(TH1D *tE, int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2){

    tE->SetMarkerStyle(marker);
    tE->SetMarkerColor(color);
    tE->SetMarkerSize(msize);
    tE->SetLineStyle(lstyle);
    tE->SetLineWidth(lwidth);
    tE->SetLineColor(color);
}

void MyPlot::StyleMe(TF1 *tE, int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2){

    tE->SetMarkerStyle(marker);
    tE->SetMarkerColor(color);
    tE->SetMarkerSize(msize);
    tE->SetLineStyle(lstyle);
    tE->SetLineWidth(lwidth);
    tE->SetLineColor(color);
}
void MyPlot::StyleMe(TGraph *tE, int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2){

    tE->SetMarkerStyle(marker);
    tE->SetMarkerColor(color);
    tE->SetMarkerSize(msize);
    tE->SetLineStyle(lstyle);
    tE->SetLineWidth(lwidth);
    tE->SetLineColor(color);
}

void MyPlot::StyleMe(TGraphErrors *tE, int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2){

	tE->SetMarkerStyle(marker);
	tE->SetMarkerColor(color);
	tE->SetMarkerSize(msize);
    tE->SetLineStyle(lstyle);
    tE->SetLineWidth(lwidth);
 	tE->SetLineColor(color);
}

// void MyPlot::StyleMe(TGraphBentErrors *tE, int marker = 20, int color = 1, float msize = 1.2, int lstyle = 1, float lwidth = 2){

// 	tE->SetMarkerStyle(marker);
// 	tE->SetMarkerColor(color);
// 	tE->SetMarkerSize(msize);
//     tE->SetLineStyle(lstyle);
//     tE->SetLineWidth(lwidth);
//  	tE->SetLineColor(color);
// }

// void MyPlot::StyleMe(TGraphAsymmErrors *tE, int marker = 20, int color = 1, float msize = 1.2, int lstyle = 1, float lwidth = 2){

// 	tE->SetMarkerStyle(marker);
// 	tE->SetMarkerColor(color);
// 	tE->SetMarkerSize(msize);
//     tE->SetLineStyle(lstyle);
//     tE->SetLineWidth(lwidth);
//  	tE->SetLineColor(color);
// }

void MyPlot::RatioBinomial(TH1D *numerator, TH1D *denominator, TH1D *ratio){

//        ratio->Divide(denominator);
//
// numerator is a subset of denominator (or vis-versa)
// all histograms must have the same binning
// error bars are calculated from binomial distribution sigma=sqrt(npq)
//
    Double_t sumN=0,sumD=0;
    Double_t p, sigma_p;        
    Double_t N, Nerr, Neff, D,Derr,Deff;

    Int_t nn = numerator->GetNbinsX();
    std::cout << "****************   " << nn << std::endl;
    for (int i=0; i<nn;i++){
      sumN=sumN+numerator->GetBinContent(i);
      sumD=sumD+denominator->GetBinContent(i); 
    }  

    if (sumN < sumD){                            // binomial distribution to get N events out of D trials, with probability N/D
      for (int i=0; i<nn; i++){
        p=0; sigma_p=0, Deff=0;
        N = numerator->GetBinContent(i);
        D = denominator->GetBinContent(i);
        Derr = denominator->GetBinError(i);      
        if (D>0) p=N/D;
        if (p>1) p=D/N;                          // safe garde needed for cases with p ~ 1
        if (Derr>0) Deff = D*D/Derr/Derr;
        if (Deff >0) sigma_p = sqrt((1-p)/p/Deff);

        if (D>0) {
            std::cout << i << "   p=" << N/D << " +/- " << N/D*sigma_p << std::endl;
            ratio->SetBinContent(i,N/D);
            ratio->SetBinError(i,N/D*sigma_p);
        }
      }
    } else {                                     // reverse ratio, binomial distribution to get D events out of N trials, with probability D/N
      for (int i=0; i<nn; i++){
        p=0; sigma_p=0, Neff=0;
        N = numerator->GetBinContent(i);
        D = denominator->GetBinContent(i);
        Nerr = numerator->GetBinError(i);      
        if (N>0) p= D/N;
        if (p>1) p= N/D;
        if (Nerr>0) Neff = N*N/Nerr/Nerr;
        if (Neff >0) sigma_p = sqrt((1-p)/p/Neff);


        if (D>0) {
            ratio->SetBinContent(i,N/D);
            ratio->SetBinError(i,N/D*sigma_p);
            std::cout << i << "   p=" << N/D << " +/- " << N/D*sigma_p << std::endl;
      }
      }
    }

     
}

