// Currently, the eta/pion ratio is smooth, having been updated.
//  Here we have eta/pion ratio for pp and eta/pion ratio for AA (RHIC level).
// In terms of the codes, pp and AA works exactly the same.

// 1. First include the file  "pp_eta_pion_ratio.C" or "AA_eta_pion_ratio.C".
// 2. Define the eta/pion ratio as an empty number array, for example y[] and the corresponding err ey[].
// 3. Define an arrary of the pt you want to have (pt[]), with the same length as that of the y and ey.
// 4. Fill pt with the numbers you want to have. And keep y[] and ey[] empty.
// 5. Run pp_eta_pion_ratio(data number, pt, y, ey)  or AA_eta_pion_ratio(...)
//    The program will fill your previous empty y and ey with the desired numbers.
#include "AA_eta_pion_ratio.C"
void long_demo()
{
  const Int_t data_num=500;

  Double_t log10pt_range[2]={-1.8,1.477};
  Double_t log10pt_step=(log10pt_range[1]-log10pt_range[0])/data_num;
  Double_t pt[data_num];
  for(Int_t i=0;i<data_num;i++)
  {
    Double_t log10pt=log10pt_range[0]+ i*log10pt_step;
    pt[i]=pow(10, log10pt );
  }
  // Double_t pt[data_num]={0.0158489,0.0160899,28.664,29.0999,29.5424};
  Double_t y[data_num];
  Double_t ey[data_num];
  Double_t ex[data_num];
  AA_eta_pion_ratio(data_num,pt,y,ey);
  for(Int_t i=0;i<data_num;i++)
  {
    ex[i]=0;
    cout<<"pt="<<pt[i]<<", "<<"y="<<y[i]<<", "<<"ey="<<ey[i]<<endl;
  }


  TCanvas* c=new TCanvas("c","c",500,1,500,500);

  c->cd();
  c->SetRightMargin(0.01);
  c->SetTopMargin(0.011);
  c->SetLogx();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.1);
  TGraphErrors* gr=new TGraphErrors(data_num,pt,y,ex,ey);
  //-------------------
  gr->GetXaxis()->SetTitleOffset(1.2);
  gr->GetXaxis()->SetTitleSize(0.04);
  gr->GetYaxis()->SetTitleOffset(1.2);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLimits(0.4,25);
  gr->SetMinimum(0);
  gr->SetMaximum(.9);
  gr->SetLineStyle(1);
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(0.1);
  gr->SetLineColor(kAzure);
  gr->SetLineWidth(2);
  gr->SetFillColorAlpha(kAzure-4,.8);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  gr->GetYaxis()->SetTitle("#eta/#pi^{0}");
  gr->Draw("E4AC");



}
