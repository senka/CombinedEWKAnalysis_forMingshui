import pyroot_logon
import limits
import os
import sys

from ROOT import *

lType = sys.argv[1]
codename = ""
planeID = sys.argv[2]

norm_sig_sm = -1
norm_bkg = -1
norm_obs = -1
if( lType == "muon" ) :
    codename = "mu"
elif( lType == "electron" ):
    codename = "el"
else:
    raise RuntimeError('InvalidLepton','You may only choose between "muon" and "electron" channels.')

basepath = '%s/src/CombinedEWKAnalysis/CommonTools/data/WV_semileptonic'%os.environ['CMSSW_BASE']

f = TFile('%s/%s_boosted.root'%(basepath,codename))

background = f.Get('background')
background_backshapeUp = f.Get('background_%sboosted_backshapeUp'%codename)
background_backshapeDown = f.Get('background_%sboosted_backshapeDown'%codename)
data_obs = f.Get('data_obs')
diboson = f.Get('diboson')

#background.Add(diboson, -1.)

norm_sig_sm = diboson.Integral()
norm_bkg = background.Integral()
norm_obs = data_obs.Integral()

theWS = RooWorkspace('WV_%sboosted'%codename, 'WV_%sboosted'%codename)

wpt = theWS.factory('W_pt[%f,%f]' % (data_obs.GetBinLowEdge(1), 
                                     data_obs.GetBinLowEdge(data_obs.GetNbinsX())+data_obs.GetBinWidth(data_obs.GetNbinsX())))
wpt.setBins(data_obs.GetNbinsX())

lz = theWS.factory('lZ[0., -1., 1.]')
# lz = theWS.factory('lZ[0.]')
lz.setConstant(False)
dkg = theWS.factory('dkg[0.,-0.15, 0.15]')
dg1 = theWS.factory('dg1[0.,-0.1,0.1]')


vars = RooArgList(wpt)
varSet = RooArgSet(wpt)

data = RooDataHist('data_obs', 'data_obs_WV_%s'%codename, vars, data_obs)
bkgHist = RooDataHist('WV_semileptonic_bkg_%s'%codename,
                      'WV_semileptonic_bkg_%s'%codename,
                      vars,
                      background)
bkgHist_systUp = RooDataHist('WV_semileptonic_bkg_%s_%sboosted_backshapeUp'%(codename,codename),
                             'WV_semileptonic_bkg_%s_%sboosted_backshapeUp'%(codename,codename),
                             vars,
                             background_backshapeUp)
bkgHist_systDown = RooDataHist('WV_semileptonic_bkg_%s_%sboosted_backshapeDown'%(codename,
                                                                                 codename),
                               'WV_semileptonic_bkg_%s_%sboosted_backshapeDown'%(codename,
                                                                                 codename),
                               vars,
                               background_backshapeDown)

dibosonHist = RooDataHist('WV_semileptonic_SM_%s_rawshape'%codename,
                          'WV_semileptonic_SM_%s_rawshape'%codename,
                          vars,
                          diboson)
dibosonPdf = RooHistFunc('WV_semileptonic_SM_%s_shape'%codename,
                         'WV_semileptonic_SM_%s_shape'%codename,
                         varSet,
                         dibosonHist)



aTGC = RooATGCFunction('ATGC_shapescale_WWgammaZ_WV_atgc_semileptonic_%s'%codename,
                       'ATGC_shapescale_%s'%codename,
                       wpt,
                       lz,
                       dkg,
                       dg1, 
                       '%s/ATGC_shape_coefficients.root'%basepath)

limtype = -1

print 'setting up for %s plane!'%planeID
if ( planeID == 'dkglZ' ):
    limtype = 0
elif ( planeID == 'dg1lZ' ):
    limtype = 1
elif ( planeID == 'dkgdg1'):
    limtype = 2
else:
    raise RuntimeError('InvalidCouplingChoice',
                       'We can only use [dkg,lZ], [dg1,lZ], and [dkg,dg1]'\
                       ' as POIs right now!')

print limtype

aTGCPdf = RooATGCSemiAnalyticPdf('ATGCPdf_WWgammaZ_WV_atgc_semileptonic_%s'%codename,
                                 'ATGCPdf_WV_%s'%codename,
                                 wpt,
                                 dkg,
                                 lz,                                 
                                 dg1,
                                 dibosonPdf,
                                 '%s/ATGC_shape_coefficients.root'%basepath,
                                 limtype)


getattr(theWS, 'import')(data)
getattr(theWS, 'import')(bkgHist)
getattr(theWS, 'import')(bkgHist_systUp)
getattr(theWS, 'import')(bkgHist_systDown)
getattr(theWS, 'import')(aTGCPdf)

theWS.Print()

fout = TFile('%s_boosted_%s_ws.root'%(codename,planeID), 'recreate')
theWS.Write()
fout.Close()

### make the card for this channel and plane ID
card = """
# Simple counting experiment, with one signal and a few background processes 
imax 1  number of channels
jmax 1  number of backgrounds
kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
shapes WV_semileptonic_bkg_{codename}  {codename}boosted ./{codename}_boosted_{planeID}_ws.root WV_{codename}boosted:$PROCESS WV_{codename}boosted:$PROCESS_$SYSTEMATIC
shapes data_obs                {codename}boosted ./{codename}_boosted_{planeID}_ws.root WV_{codename}boosted:$PROCESS
shapes WWgammaZ_WV_atgc_semileptonic_{codename} {codename}boosted ./{codename}_boosted_{planeID}_ws.root WV_{codename}boosted:ATGCPdf_$PROCESS
------------
bin {codename}boosted 
observation {norm_obs}
------------
bin                         {codename}boosted		      {codename}boosted
process                     WWgammaZ_WV_atgc_semileptonic_{codename}   WV_semileptonic_bkg_{codename}
process                     0			      1		
rate                        {norm_sig_sm}		      {norm_bkg}	

------------
lumi_8TeV           lnN     1.044		      -
CMS_eff_{codename[0]}           lnN     1.02                      -
CMS_trigger_{codename[0]}       lnN     1.01                      -
{codename}boosted_backshape shape1  -                         1.0 
sigXSsyst           lnN     1.034                     -
""".format(codename=codename,planeID=planeID,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)

print card

cardfile = open('wv_semil_%sboosted_%s.txt'%(codename,planeID),'w')
cardfile.write(card)
cardfile.close
