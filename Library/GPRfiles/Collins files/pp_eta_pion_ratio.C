// Yuanjie Ren August 2, 2020
const double x_range[2]={-1.8,1.477};
const int nPred=500;
Int_t pp_eta_pion_ratio(const Int_t data_number, Double_t input_pt[], Double_t input_EtaPionRatio[], Double_t input_EtaPionRatioErr[])
{
  ifstream filein;
  string trashline;
  int output=0;


  double pp_pt[nPred],pp_y[nPred],pp_ey[nPred];
  filein.open("pp_maximal_range_smooth.txt");
  int header_lines=2;
  for(int i=0;i<header_lines;i++)
  {
    getline(filein,trashline);
  if(output==1)  cout<<"trash: "<<trashline<<endl;
  }
  for(int i=0;i<nPred;i++)
  {
    filein>>pp_pt[i]>>pp_y[i]>>pp_ey[i];
    if(output==1)cout<<"i="<<i<<", (x,y,ey)=("<<pp_pt[i]<<","<<pp_y[i]<<","<<pp_ey[i]<<")"<<endl;
  }
  filein.close();



  for(Int_t i=0;i<data_number;i++)
  {
    input_EtaPionRatio[i]=continuity(input_pt[i],pp_y);
    input_EtaPionRatioErr[i]=continuity(input_pt[i],pp_ey);
    // cout<<"y="<<input_EtaPionRatio[i]<<endl;
    // cout<<"ey="<<input_EtaPionRatioErr[i]<<endl;
  }




return 0;

}//end of main function


      double continuity(double pt, double quantity[])
       {
         double precise_index=Target_index(pt);
         int index=floor(Target_index( pt));
         double alpha=precise_index-index ;
         return quantity[index]*(1-alpha)+quantity[index+1]*alpha;
       }

      //----------------------------------

      double Target_index(double pt)
      {
        double x_upper_limit=x_range[1];
        double x_lower_limit=x_range[0];
        double width=x_upper_limit-x_lower_limit;
        double gap=width/nPred;
        return (log10(pt)-x_lower_limit)/gap;
      }
