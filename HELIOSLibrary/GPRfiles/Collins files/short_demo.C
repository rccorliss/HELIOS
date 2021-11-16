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
void short_demo()
{
  const Int_t data_num=6;
  Double_t pt[data_num]={0.0158489,0.0160899,0.0160899,28.664,29.0999,29.5424};
  Double_t y[data_num];
  Double_t ey[data_num];
  AA_eta_pion_ratio(data_num,pt,y,ey);
  for(Int_t i=0;i<data_num;i++)
  {
    cout<<"pt="<<pt[i]<<", "<<"y="<<y[i]<<", "<<"ey="<<ey[i]<<endl;
  }




}
