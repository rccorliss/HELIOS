//
// Template functions to make nice looking plot
//
// Axel Drees 1/1/2015
//            9/17/2021 converted to class structure
//
// header files contains all definitions and default values for papameters
// 
//
#ifndef MyPlot_h
#define MyPlot_h

class MyPlot {

public:
MyPlot();                                                // default constructor
virtual ~MyPlot();        

void RatioBinomial(TH1D *numerator, TH1D *denominator, TH1D *ratio);

TCanvas *Canvas(const Char_t *name, Int_t px2=720, Int_t py2 =920, Int_t px1=10,Int_t py1=10, Int_t logy=0, Int_t logx=0 ); 
TH1D    *Frame(const Char_t *name, const Char_t *xaxis, const Char_t *yaxis, Double_t xmin=0., Double_t xmax=10., Double_t ymin=0., Double_t ymax=10., Int_t center=0);
TLegend *Legend(const Char_t *name, Double_t x1=0.6, Double_t y1=0.79,Double_t x2=0.8,Double_t y2=0.94, Double_t size=1.);

void StyleMe (TH1D *tE,              int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2);
void StyleMe (TF1 *tE,               int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2);
void StyleMe (TGraph *tE,            int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2);
void StyleMe (TGraphErrors *tE,      int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2);
//void StyleMe (TGraphBentErrors *tE,  int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2);
//void StyleMe (TGraphAsymmErrors *tE, int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2);

void SetFont(int i) {Font = i;}
void SetxTitleOffset(float t) {xTitleOffset = t;}
void SetyTitleOffset(float t) {yTitleOffset = t;}
void SetTitleSize(float t) {TitleSize = t;}
void SetLabelSize(float t) {LabelSize = t;}
void SetLegendSize(float t) {LegendSize = t;}
void SetLeftMargin(float t) {LeftMargin = t;}
void SetRightMargin(float t) {RightMargin = t;}
void SetBottomMargin(float t) {BottomMargin = t;}
void SetTopMargin(float t) {TopMargin = t;}
void SetTitleColor(int color) {TitleColor = color;}
void SetLegendColor(int color) {LegendColor = color;}

void Reset(){
    Font = FontDefault;
    TitleSize = TitleSizeDefault;
    xTitleOffset = xTitleOffsetDefault;
    yTitleOffset = yTitleOffsetDefault;
    LabelSize = LabelSizeDefault;
    LegendSize = LegendSizeDefault;
    LeftMargin = LeftMarginDefault;
    RightMargin = RightMarginDefault;
    BottomMargin = BottomMarginDefault;
    TopMargin = TopMarginDefault;  
    TitleColor = kBlack;
    LegendColor = kBlack;
  }


private:


const int FontDefault = 42;
const float TitleSizeDefault = 0.06;
const float xTitleOffsetDefault = 1.1;
const float yTitleOffsetDefault = 1.1;
const float LabelSizeDefault = 0.05;
const float LegendSizeDefault = 0.04;
const float LeftMarginDefault = 0.18;
const float RightMarginDefault = 0.02;
const float BottomMarginDefault = 0.15;
const float TopMarginDefault = 0.02;

int Font = FontDefault;
float TitleSize = TitleSizeDefault;
float xTitleOffset = xTitleOffsetDefault;
float yTitleOffset = yTitleOffsetDefault;
float LabelSize = LabelSizeDefault;
float LegendSize = LegendSizeDefault;
float LeftMargin = LeftMarginDefault;
float RightMargin = RightMarginDefault;
float BottomMargin = BottomMarginDefault;
float TopMargin = TopMarginDefault;  
float   TitleColor = kBlack;
float   LegendColor = kBlack;

}; // end of Particle class

#endif

