import pyroot_logon
import limits

from ROOT import *

basepath = '/home/lgray/CMSSW_6_2_0_pre5/src/CombinedEWKAnalysis/CommonTools/data/WV_semileptonic'

f = TFile('%s/mu_boosted.root'%basepath)

background = f.Get('background')
data_obs = f.Get('data_obs')
diboson = f.Get('diboson')

background.Add(diboson, -1.)

theWS = RooWorkspace('theWS', 'theWS')

wpt = theWS.factory('W_pt[%f,%f]' % (data_obs.GetBinLowEdge(1), 
                                     data_obs.GetBinLowEdge(data_obs.GetNbinsX())+data_obs.GetBinWidth(data_obs.GetNbinsX())))
wpt.setBins(data_obs.GetNbinsX())

lz = theWS.factory('lZ[0., -0.03, 0.03]')
# lz = theWS.factory('lZ[0.]')
lz.setConstant(False)
dkg = theWS.factory('dkg[0., -0.15, 0.15]')
dg1 = theWS.factory('dg1[0.]')


vars = RooArgList(wpt)
varSet = RooArgSet(wpt)

data = RooDataHist('data_obs', 'data_obs', vars, data_obs)
bkgHist = RooDataHist('WV_semileptonic_bkg_mu',
                      'WV_semileptonic_bkg_mu',
                      vars, background)
dibosonHist = RooDataHist('WV_atgc_semileptonic_mu',
                          'WV_atgc_semileptonic_mu',
                          vars, diboson)

bkgPdf = RooHistPdf('bkgPdf', 'bkgPdf', varSet, bkgHist)
dibosonPdf = RooHistPdf('dibosonPdf', 'dibosonPdf', varSet, dibosonHist,0)
dibosonFunc = RooHistFunc('dibosonFunc', 'dibosonPdf', varSet, dibosonHist,0)

print '%s/ATGC_shape_coefficients.root'%basepath

aTGC = RooATGCFunction('aTGC', 'aTGC', wpt, lz, dkg, dg1, 
                       '%s/ATGC_shape_coefficients.root'%basepath)

aTGCPdfOld = RooGenericPdf('aTGCPdfOld','aTGCPdfOld','@0*@1',RooArgList(dibosonPdf,aTGC))

aTGCPdf = RooATGCSemiAnalyticPdf('aTGCPdf',
                                 'aTGCPdf_WV_mu',
                                 wpt,
                                 dkg,
                                 lz,                                 
                                 dg1,
                                 dibosonFunc,
                                 '%s/ATGC_shape_coefficients.root'%basepath)

nbkg = theWS.factory('prod::bkg_yield(n_bkg[%f],bkg_nrm[1.,-1,5])' % \
                         (background.Integral()))
ndiboson = theWS.factory('prod::diboson_yield(n_diboson[%f],diboson_nrm[1.,-1,5])' % \
                             (diboson.Integral()))


print "starting normalizations"
print nbkg.getVal()
print ndiboson.getVal()

getattr(theWS, 'import')(data)
getattr(theWS, 'import')(bkgPdf)
getattr(theWS, 'import')(aTGCPdfOld)
getattr(theWS, 'import')(aTGCPdf, RooFit.RecycleConflictNodes())

theWS.factory('RooExtendPdf::bkg_extended(bkgPdf,bkg_yield)')
theWS.factory('RooExtendPdf::aTGC_extended(aTGCPdf,diboson_yield)')
comps = RooArgList(theWS.argSet('aTGC_extended,bkg_extended'))
total = RooAddPdf('total', 'total', comps)

getattr(theWS, 'import')(total)
theWS.factory('RooGaussian::bkg_const(bkg_nrm, 1.0, 0.05)')
total_const = theWS.factory('PROD::total_const(total, bkg_const)')

theWS.factory('RooATGCProcessScaling::testScale(W_pt,dkg,lZ,dg1,dibosonPdf,"%s/ATGC_shape_coefficients.root")'%basepath)



theWS.Print()

total_const.fitTo(data, RooFit.Extended())

frame = wpt.frame()
data.plotOn(frame)
total.plotOn(frame)
total.plotOn(frame, RooFit.Components('bkg*'),
                   RooFit.LineColor(kRed),
                   RooFit.LineStyle(kDashed))
data.plotOn(frame)
frame.Draw()

gPad.Update()
#gPad.WaitPrimitive()

poi = RooArgSet(lz, dkg)

theWS.var('dkg').setVal(0.0)
theWS.var('lZ').setVal(0.0)

print theWS.function('testScale').getVal()

theWS.var('dkg').setVal(0.1)

print theWS.function('testScale').getVal()

#limit = limits.plcLimit(wpt, poi, total_const, theWS, data, verbose = True)

#theWS.Print()
#fout = TFile('ATGC_likelihood.root', 'recreate')
#theWS.Write()
#fout.Close()
