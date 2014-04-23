{
gROOT->Reset();

gSystem->Load("libFWCoreFWLite.so");
AutoLibraryLoader::enable();
gSystem->Load("libDataFormatsFWLite.so");

gSystem->Load("${CMSSW_BASE}/lib/${SCRAM_ARCH}/libHiggsAnalysisCombinedLimit.so");


//Suppress the printout from making plots
 gErrorIgnoreLevel = 2000;

RooFit::RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling) ;

  TString path = gSystem->GetIncludePath();
  path += "-I. -I$ROOTSYS/src -I$ROOFITSYS/include";
  gSystem->SetIncludePath(path.Data());


}
