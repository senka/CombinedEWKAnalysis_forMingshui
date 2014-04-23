void dealWith_unEqualBinning_ZZ(TString file_name){

  TFile * file_in=new TFile(file_name,"read");
  //  file_in->cd();
  file_in->cd("ll");
  //  ll->cd();
  file_in.ls();
  cout << "listing dir.."<<endl;
  gDirectory->ls();
  gDirectory->cd();
  //  TList * input_p=file_in->GetList();
  //  TList * list_histos=file_in->GetListOfKeys();
  TList * list_histos=gDirectory->GetListOfKeys();
  cout << "--------> histo name: "<< list_histos->At(0)->GetName()<< endl;
  TString file_out=file_name.Remove(file_name.Length() - 5, 5);
  file_out+="_uniformBins.root";
  cout <<" ----> output: "<< file_out<<"\n" << endl;
  TFile *outFile=new TFile(file_out,"recreate");
  outFile->cd();

  for (int i=0;i<list_histos->GetEntries();i++){
    file_in->cd("ll");
    TString histo_name= list_histos->At(i)->GetName();

    //    TH1D * plot=file_in->Get(histo_name);
    TH1D * plot=gDirectory->Get(histo_name);
    //cout << " plot \""<<histo_name<<"\""<< endl;
    cout <<"\treading histo "<< histo_name<< " and rebinning.." << endl;
    cout <<"\t\t integral: "<<plot->Integral() << endl;
    int N_bins=plot->GetNbinsX();
    TH1D * plot_rebin=new TH1D("plot_rebin","plot_rebin",N_bins-2,0,N_bins-2);
    for (int j=1;j<=N_bins;j++)
      plot_rebin->SetBinContent(j,plot->GetBinContent(j+2));
    cout <<"\t\tbin yield starting: "<< plot->Integral()<<" after rebinning: "<<plot_rebin->Integral()  << endl;

    outFile->cd();
    plot_rebin->Write(histo_name);
    delete plot_rebin;
  }

  outFile->Close();
  file_in->Close();

  /*
  list_histos->GetEntries();
  list_histos->At(4)->GetName();
  TH1D * plot=file_in->Get("data_obs");
  */

}
