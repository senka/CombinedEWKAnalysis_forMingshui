import pyroot_logon
import limits
import os
import sys

from array import *


from ROOT import *
from optparse import OptionParser
from ConfigParser import SafeConfigParser


def isItCorrelated(name):
    print '\t ----> isItCorrelated: testing ',name
    if ('_eff_b' in name or '_les' in name or '_pu' in name or '_umet' in name or '_res_j' in name or '_scale_j' in name ):
        print '-> true'
        return True
    else:
        print '-> false'
        return False


parser = OptionParser(description="%prog : A RooStats Implementation of Anomalous Triple Gauge Coupling Analysis.",
                      usage="buildWZworkspace --config=example_config.cfg")
cfgparse = SafeConfigParser()

parser.add_option("--config",dest="config",help="The name of the input configuration file.")
(options,args) = parser.parse_args()

miss_options = False

if options.config is None:
    print 'Need to specify --config'
    miss_options=True
    
if miss_options:
    exit(1)
        
cfgparse.read(options.config)
options.config = cfgparse # put the parsed config file into our options

cfg = options.config


#lType = sys.argv[1]
#codename = ""
#planeID = sys.argv[2]

norm_sig_sm = -1
norm_sig_sm_up = -1
norm_sig_sm_down = -1
norm_bkg = -1
norm_obs = -1


fit_sections = cfg.sections()
fit_sections.remove('Global') #don't need to iterate over the global configuration

basepath = '%s/src/CombinedEWKAnalysis/CommonTools/data/WV_semileptonic'%os.environ['CMSSW_BASE']

for section in fit_sections:
    codename = section
    lType = codename
    print '\n\tlType=',lType
    f = TFile('%s/%s_boosted.root'%(basepath,codename))

    Nbkg = cfg.get(codename,'Nbkg')
    print "Nbkg= ",Nbkg
    Nbkg_int=int(Nbkg)

    bkg_name = []
    for i in range(1,Nbkg_int+1):
        bkg_name.append(cfg.get(codename,'bkg%i_name'%i))

    background = []

    for i in range(0,Nbkg_int):
        background.append(f.Get(bkg_name[i]))


    print 'backgrounds= ',background
    background_shapeSyst = []
    for i in range(0,Nbkg_int):
        background_shapeSyst.append([])
        for name in cfg.get(codename,'bkg%i_shape_syst'%(i+1)).split(','):
            background_shapeSyst[i].append(name)



    background_backshapeUp = []
    background_backshapeDown = []

    for j in range(0,Nbkg_int):
        background_backshapeUp.append([])
        background_backshapeDown.append([])
        for i in range(0,len(background_shapeSyst[j])):
            print ' bkg shape syst: ',background_shapeSyst[j]
            print ' getting bkgUp ','%sUp'%background_shapeSyst[j][i]
            background_backshapeUp[j].append(f.Get('%sUp'%background_shapeSyst[j][i]))
            background_backshapeDown[j].append(f.Get('%sDown'%background_shapeSyst[j][i]))


    data_obs = f.Get('data_obs')
#    diboson = f.Get('diboson')
    diboson = f.Get('zz2l2nu')

    doSignalShape_unc=False
    cfg_items=cfg.items(codename)
    for cfg_item in cfg_items:
        if 'signal_shape_syst' in cfg_item:
            doSignalShape_unc = True

    print 'doSignalShape_unc=',doSignalShape_unc

    if (doSignalShape_unc):
        diboson_up = {}
        diboson_down = {}
        norm_sig_sm_up = {}
        norm_sig_sm_down = {}
        signal_shapeSyst = [string(i) for i in cfg.get(codename,'signal_shape_syst').split(',')]
        for i in range(0,len(signal_shapeSyst)):
            print ' signal shape syst: ',signal_shapeSyst[i]
            diboson_up[i] = f.Get('%sUp'%signal_shapeSyst[i])
            diboson_down[i] = f.Get('%sDown'%signal_shapeSyst[i])
            norm_sig_sm_up[i] = diboson_up[i].Integral()
            norm_sig_sm_down[i] = diboson_down[i].Integral()


    
    norm_sig_sm = diboson.Integral()
    norm_bkg = []
    for i in range(0,Nbkg_int):
        norm_bkg.append(background[i].Integral())
    norm_obs = data_obs.Integral()
    print 'bkg integral: ',norm_bkg
    
    if (doSignalShape_unc):
        print 'signal shape unc: ',norm_sig_sm_down,' ',norm_sig_sm,' ',norm_sig_sm_up
    
    theWS = RooWorkspace('WV_%sboosted'%codename, 'WV_%sboosted'%codename)
    
    wpt = theWS.factory('W_pt_%s[%f,%f]' % (codename,data_obs.GetBinLowEdge(1), 
                                         data_obs.GetBinLowEdge(data_obs.GetNbinsX())+data_obs.GetBinWidth(data_obs.GetNbinsX())))

    binning=array('d',[])

    for i in range(1, data_obs.GetNbinsX()+1):
        binning.append(data_obs.GetBinLowEdge(i))
    binning.append(data_obs.GetBinLowEdge(data_obs.GetNbinsX()+1))

    print "bining: "
    for i in range(0, len(binning)):
        print binning[i]

    bins=RooBinning(len(binning)-1, binning)

    wpt.setBinning(bins)
    lz = theWS.factory('lZ[0., -0.006, 0.006]')
    lz.setConstant(False)
    dkg = theWS.factory('dkg[0.,-0.006, 0.006]')
    dg1 = theWS.factory('dg1[0.,-0.006,0.006]')


    vars = RooArgList(wpt)
    varSet = RooArgSet(wpt)

    data = RooDataHist('data_obs', 'data_obs_WV_%s'%codename, vars, data_obs)
    bkgHist = {}
    for i in range(0,Nbkg_int):
        bkgHist[i] = RooDataHist('WV_semileptonic_bkg%i_%s'%(i+1,codename),
                              'WV_semileptonic_bkg%i_%s'%(i+1,codename),
                              vars,
                              background[i])
    
    bkgHist_systUp = []
    bkgHist_systDown = []

    for j in range(0,Nbkg_int):
        bkgHist_systUp.append([])
        bkgHist_systDown.append([])
        for i in range(0,len(background_shapeSyst[j])):
            print j," ",i
            #change name here:
            print '\n\t\t==========> testing in function= ', isItCorrelated('testing')
            print '\t\t==========> wz3lnu_CMS_eff_b in function= ', isItCorrelated('wz3lnu_CMS_eff_b')
            print '\t\t==========> %s in function= '%background_shapeSyst[j][i], isItCorrelated(background_shapeSyst[j][i])
            print '\n'
            if (isItCorrelated(background_shapeSyst[j][i])):
                print ' \n\t\t ==================================> <=========================== '
                name_forCorr=background_shapeSyst[j][i]
                print ' name_forCorr= ',name_forCorr
                name_forCorr=name_forCorr.replace('zll_','')
                name_forCorr=name_forCorr.replace('wz3lnu_','')
                print ' -> name_forCorr= ',name_forCorr
                bkgHist_systUp[j].append(RooDataHist('WV_semileptonic_bkg%i_%s_%sUp'%(j+1,codename,name_forCorr),
                                                    'WV_semileptonic_bkg%i_%s_%sUp'%(j+1,codename,name_forCorr),
                                                    vars,
                                                    background_backshapeUp[j][i]))
                bkgHist_systDown[j].append(RooDataHist('WV_semileptonic_bkg%i_%s_%sDown'%(j+1,codename,name_forCorr),
                                                       'WV_semileptonic_bkg%i_%s_%sDown'%(j+1,codename,name_forCorr),
                                                       vars,
                                                       background_backshapeDown[j][i]))
            else:
                bkgHist_systUp[j].append(RooDataHist('WV_semileptonic_bkg%i_%s_%sUp'%(j+1,codename,background_shapeSyst[j][i]),
                                                     'WV_semileptonic_bkg%i_%s_%sUp'%(j+1,codename,background_shapeSyst[j][i]),
                                                     vars,
                                                     background_backshapeUp[j][i]))
                bkgHist_systDown[j].append(RooDataHist('WV_semileptonic_bkg%i_%s_%sDown'%(j+1,codename,background_shapeSyst[j][i]),
                                                       'WV_semileptonic_bkg%i_%s_%sDown'%(j+1,codename,background_shapeSyst[j][i]),
                                                       vars,
                                                       background_backshapeDown[j][i]))

#            bkgHist_systUp[j].append(RooDataHist('WV_semileptonic_bkg%i_%s_%sUp'%(j+1,codename,background_shapeSyst[j][i]),
#                                                 'WV_semileptonic_bkg%i_%s_%sUp'%(j+1,codename,background_shapeSyst[j][i]),
#                                                 vars,
#                                                 background_backshapeUp[j][i]))
#            bkgHist_systDown[j].append(RooDataHist('WV_semileptonic_bkg%i_%s_%sDown'%(j+1,codename,background_shapeSyst[j][i]),
#                                                   'WV_semileptonic_bkg%i_%s_%sDown'%(j+1,codename,background_shapeSyst[j][i]),
#                                                   vars,
#                                                   background_backshapeDown[j][i]))

    dibosonHist = RooDataHist('WV_semileptonic_SM_%s_rawshape'%codename,
                              'WV_semileptonic_SM_%s_rawshape'%codename,
                              vars,
                              diboson)

    if (doSignalShape_unc):
        dibosonHist_up = {}
        dibosonHist_down = {}
        for i in range(0,len(signal_shapeSyst)):
            print ' \n\t\t ==================================> SIGNAL  %s <=========================== '%signal_shapeSyst[i],isItCorrelated(str(signal_shapeSyst[i]))
            print ' \n\t\t ==================================> SIGNAL  zz2l2nu_CMS_scale_j <=========================== ',isItCorrelated('zz2l2nu_CMS_scale_j')
            #change name here
#            dibosonHist_up[i] = RooDataHist('WV_semileptonic_SM_%s_rawshape_%sUp'%(codename,signal_shapeSyst[i]),
 #                                        'WV_semileptonic_SM_%s_rawshape_%sUp'%(codename,signal_shapeSyst[i]),
  #                                       vars,
   #                                      diboson_up[i])
    #        dibosonHist_down[i] = RooDataHist('WV_semileptonic_SM_%s_rawshape_%sDown'%(codename,signal_shapeSyst[i]),
     #                                      'WV_semileptonic_SM_%s_rawshape_%sDown'%(codename,signal_shapeSyst[i]),
      #                                     vars,
       #                                    diboson_down[i])

            if (isItCorrelated(str(signal_shapeSyst[i]))):
                print ' \n\t\t ==================================> <=========================== '
                name_forCorr=str(signal_shapeSyst[i])
                print ' name_forCorr= ',name_forCorr
                name_forCorr=name_forCorr.replace('zz2l2nu_','')
                print ' -> name_forCorr= ',name_forCorr
                dibosonHist_up[i] = RooDataHist('WV_semileptonic_SM_%s_rawshape_%sUp'%(codename,name_forCorr),
                                                'WV_semileptonic_SM_%s_rawshape_%sUp'%(codename,name_forCorr),
                                                vars,
                                                diboson_up[i])
                dibosonHist_down[i] = RooDataHist('WV_semileptonic_SM_%s_rawshape_%sDown'%(codename,name_forCorr),
                                                  'WV_semileptonic_SM_%s_rawshape_%sDown'%(codename,name_forCorr),
                                                  vars,
                                                  diboson_down[i])
            else:
                dibosonHist_up[i] = RooDataHist('WV_semileptonic_SM_%s_rawshape_%sUp'%(codename,signal_shapeSyst[i]),
                                                'WV_semileptonic_SM_%s_rawshape_%sUp'%(codename,signal_shapeSyst[i]),
                                                vars,
                                                diboson_up[i])
                dibosonHist_down[i] = RooDataHist('WV_semileptonic_SM_%s_rawshape_%sDown'%(codename,signal_shapeSyst[i]),
                                                  'WV_semileptonic_SM_%s_rawshape_%sDown'%(codename,signal_shapeSyst[i]),
                                                  vars,
                                                  diboson_down[i])
                
            

    dibosonPdf = RooHistFunc('WV_semileptonic_SM_%s_shape'%codename,
                             'WV_semileptonic_SM_%s_shape'%codename,
                             varSet,
                             dibosonHist)

    if (doSignalShape_unc):
        dibosonPdf_up = {}
        dibosonPdf_down = {}
        for i in range(0,len(signal_shapeSyst)):
            # change name here
            if (isItCorrelated(str(signal_shapeSyst[i]))):
                name_forCorr=str(signal_shapeSyst[i])
                name_forCorr=name_forCorr.replace('zz2l2nu_','')
                dibosonPdf_up[i] = RooHistFunc('WV_semileptonic_SM_%s_shape_%sUp'%(codename,name_forCorr),
                                               'WV_semileptonic_SM_%s_shape_%sUp'%(codename,name_forCorr),
                                               varSet,
                                               dibosonHist_up[i])
                dibosonPdf_down[i] = RooHistFunc('WV_semileptonic_SM_%s_shape_%sDown'%(codename,name_forCorr),
                                                 'WV_semileptonic_SM_%s_shape_%sDown'%(codename,name_forCorr),
                                                 varSet,
                                                 dibosonHist_down[i])

            else:
                dibosonPdf_up[i] = RooHistFunc('WV_semileptonic_SM_%s_shape_%sUp'%(codename,signal_shapeSyst[i]),
                                               'WV_semileptonic_SM_%s_shape_%sUp'%(codename,signal_shapeSyst[i]),
                                               varSet,
                                               dibosonHist_up[i])
                dibosonPdf_down[i] = RooHistFunc('WV_semileptonic_SM_%s_shape_%sDown'%(codename,signal_shapeSyst[i]),
                                                 'WV_semileptonic_SM_%s_shape_%sDown'%(codename,signal_shapeSyst[i]),
                                                 varSet,
                                                 dibosonHist_down[i])
    

#    print '\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ reading RooATGCFunction\n'
#    aTGC = RooATGCFunction_wz('ATGC_shapescale_WWgammaZ_WV_atgc_semileptonic_%s'%codename,
#                              'ATGC_shapescale_%s'%codename,
#                              wpt,
#                              lz,
#                              dkg,
#                              dg1, 
#                              '%s/signal_WV_%s.root'%(basepath,codename))
    
    print '\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ read RooATGCFunction\n'
    limtype = -1

    
    planeID = 'dkglZ'
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

    print '\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ reading RooATGCSemi\n'

    if (doSignalShape_unc):
        kappaLow = {}
        kappaHigh = {}
        aTGCPdf_norm = {}
        theta = {}
        kappaLow_sum_d = 1.
        kappaHigh_sum_d = 1.
        
        for i in range(0,len(signal_shapeSyst)):
            kappaLow[i] = RooRealVar("kappaL_%s"%(i+1),"kappaL_%s"%(i+1),norm_sig_sm_down[i]/norm_sig_sm)
            kappaLow[i].setConstant(True)
            kappaHigh[i] = RooRealVar("kappaH_%s"%(i+1),"kappaH_%s"%(i+1),norm_sig_sm_up[i]/norm_sig_sm)
            kappaHigh[i].setConstant(True)
            kappaLow_sum_d = kappaLow_sum_d*norm_sig_sm_down[i]/norm_sig_sm
            kappaHigh_sum_d = kappaHigh_sum_d*norm_sig_sm_up[i]/norm_sig_sm
            
#            theWS.factory("%s[-7,7]"%signal_shapeSyst[i])
#            theta[i] = theWS.var("%s"%signal_shapeSyst[i])

            if (isItCorrelated(str(signal_shapeSyst[i]))):
                name_forCorr=str(signal_shapeSyst[i])
                name_forCorr=name_forCorr.replace('zz2l2nu_','')
                theWS.factory("%s[-7,7]"%name_forCorr)
                theta[i] = theWS.var("%s"%name_forCorr)

            else:
                theWS.factory("%s[-7,7]"%signal_shapeSyst[i])
                theta[i] = theWS.var("%s"%signal_shapeSyst[i])
            
            aTGCPdf_norm[i] = AsymPow('ATGCPdf_WWgammaZ_WV_atgc_semileptonic_%s_integral%s'%(codename,i+1),
                                      'ATGCPdf_WV_%s_integral%s'%(codename,i+1),
                                      kappaLow[i],
                                      kappaHigh[i],
                                      theta[i])

        if (len(signal_shapeSyst)==1):
            aTGCPdf_norm_sum = aTGCPdf_norm[0]
        else:
            for i in range(0,len(signal_shapeSyst)):
                if (i==0): prodset=RooArgList(aTGCPdf_norm[i])
                else: prodset.add(RooArgList(aTGCPdf_norm[i]))
            aTGCPdf_norm_sum = RooProduct("aTGCPdf_norm_sum","aTGCPdf_norm_sum",prodset)

        kappaLow_sum = RooRealVar("kappaLow_sum","kappaLow_sum",kappaLow_sum_d)
        kappaHigh_sum = RooRealVar("kappaHigh_sum","kappaHigh_sum",kappaHigh_sum_d)

        aTGCPdf_norm_sum.SetNameTitle('ATGCPdf_WWgammaZ_WV_atgc_semileptonic_%s_norm'%codename,
                                      'ATGCPdf_WV_%s_norm'%codename)
        
        
    aTGCPdf = RooATGCSemiAnalyticPdf_wz('ATGCPdf_WWgammaZ_WV_atgc_semileptonic_%s'%codename,
                                        'ATGCPdf_WV_%s'%codename,
                                        wpt,
                                        dkg,
                                        lz,                                 
                                        dg1,
                                        dibosonPdf,
                                        '%s/signal_WV_%s.root'%(basepath,codename),
                                        limtype
                                        )

    if (doSignalShape_unc):
        aTGCPdf_up = {}
        aTGCPdf_down = {}
        for i in range(0,len(signal_shapeSyst)):
            # change name here
            if (isItCorrelated(str(signal_shapeSyst[i]))):
                name_forCorr=str(signal_shapeSyst[i])
                name_forCorr=name_forCorr.replace('zz2l2nu_','')
                aTGCPdf_up[i] = RooATGCSemiAnalyticPdf_wz('ATGCPdf_WWgammaZ_WV_atgc_semileptonic_%s_%sUp'%(codename,name_forCorr),
                                                          'ATGCPdf_WV_%s'%codename,
                                                          wpt,
                                                          dkg,
                                                          lz,                                 
                                                          dg1,
                                                          dibosonPdf_up[i],
                                                          '%s/signal_WV_%s.root'%(basepath,codename),
                                                          limtype
                                                          )
                aTGCPdf_down[i] = RooATGCSemiAnalyticPdf_wz('ATGCPdf_WWgammaZ_WV_atgc_semileptonic_%s_%sDown'%(codename,name_forCorr),
                                                            'ATGCPdf_WV_%s'%codename,
                                                            wpt,
                                                            dkg,
                                                            lz,                                 
                                                            dg1,
                                                            dibosonPdf_down[i],
                                                            '%s/signal_WV_%s.root'%(basepath,codename),
                                                            limtype
                                                            )

            else:
                aTGCPdf_up[i] = RooATGCSemiAnalyticPdf_wz('ATGCPdf_WWgammaZ_WV_atgc_semileptonic_%s_%sUp'%(codename,signal_shapeSyst[i]),
                                                          'ATGCPdf_WV_%s'%codename,
                                                          wpt,
                                                          dkg,
                                                          lz,                                 
                                                          dg1,
                                                          dibosonPdf_up[i],
                                                          '%s/signal_WV_%s.root'%(basepath,codename),
                                                          limtype
                                                          )
                aTGCPdf_down[i] = RooATGCSemiAnalyticPdf_wz('ATGCPdf_WWgammaZ_WV_atgc_semileptonic_%s_%sDown'%(codename,signal_shapeSyst[i]),
                                                            'ATGCPdf_WV_%s'%codename,
                                                            wpt,
                                                            dkg,
                                                            lz,                                 
                                                            dg1,
                                                            dibosonPdf_down[i],
                                                            '%s/signal_WV_%s.root'%(basepath,codename),
                                                            limtype
                                                            )
    
    
    print '\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ read RooATGCSemi\n'
    
    getattr(theWS, 'import')(data)
    for i in range(0,Nbkg_int):
        getattr(theWS, 'import')(bkgHist[i])
    for j in range(0,Nbkg_int):
        for i in range(0,len(background_shapeSyst[j])):
            getattr(theWS, 'import')(bkgHist_systUp[j][i])
            getattr(theWS, 'import')(bkgHist_systDown[j][i])
    getattr(theWS, 'import')(aTGCPdf)
    if (doSignalShape_unc):
        for i in range(0,len(signal_shapeSyst)):
            getattr(theWS, 'import')(aTGCPdf_up[i])
            getattr(theWS, 'import')(aTGCPdf_down[i])
#            getattr(theWS, 'import')(aTGCPdf_norm[i])
        getattr(theWS, 'import')(aTGCPdf_norm_sum)
    
    theWS.Print()
    
    fout = TFile('%s_boosted_ws.root'%(codename), 'recreate')
    theWS.Write()
    fout.Close()

### make the card for this channel and plane ID
    card = """
# Simple counting experiment, with one signal and a few background processes 
imax 1  number of channels
jmax {Nbkg_int}  number of backgrounds
kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,Nbkg_int=Nbkg_int)
    for i in range(0,Nbkg_int):
        card += """
shapes WV_semileptonic_bkg{Nbkg_int}_{codename}  {codename}boosted ./{codename}_boosted_ws.root WV_{codename}boosted:$PROCESS WV_{codename}boosted:$PROCESS_$SYSTEMATIC""".format(Nbkg_int=i+1,codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
    card += """
shapes data_obs                {codename}boosted ./{codename}_boosted_ws.root WV_{codename}boosted:$PROCESS """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,Nbkg_int=Nbkg_int)
    if (doSignalShape_unc):
        card += """   
shapes WWgammaZ_WV_atgc_semileptonic_{codename} {codename}boosted ./{codename}_boosted_ws.root WV_{codename}boosted:ATGCPdf_$PROCESS WV_{codename}boosted:ATGCPdf_$PROCESS_$SYSTEMATIC """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
    else:
        card += """   
shapes WWgammaZ_WV_atgc_semileptonic_{codename} {codename}boosted ./{codename}_boosted_ws.root WV_{codename}boosted:ATGCPdf_$PROCESS
""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
        
    card += """   
------------
bin {codename}boosted 
observation {norm_obs}
------------
bin                         {codename}boosted\t\t""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)

    for i in range(0,Nbkg_int):
        card += """\t\t\t{codename}boosted""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    card += """       
process\t\t\t    WWgammaZ_WV_atgc_semileptonic_{codename}    """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    for i in range(0,Nbkg_int):
        card += """\tWV_semileptonic_bkg{Nbkg_int}_{codename}""".format(Nbkg_int=i+1,codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    card += """       
process                     0	  	      """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    for i in range(0,Nbkg_int):
        card += """ \t\t\t\t{i}""".format(i=i+1,codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)
        
    card += """       
rate                        {norm_sig_sm}\t""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    for i in range(0,Nbkg_int):
        card += """ \t\t\t{norm_bkg}""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

	
    card += """           
------------
lumi_8TeV       \t    lnN \t    1.026             """
    for i in range(0,Nbkg_int):
        if (i==0 or i==2): card += """\t\t\t\t1.026""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
        else: card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
    
    card += """       
CMS_eff_{codename[0]}           lnN     1.03          """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg[i],norm_obs=norm_obs)

    for i in range(0,Nbkg_int):
        if (i==0 or i==2): card += """\t\t\t\t1.03""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
        else: card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)

    card += """       
CMS_hzz2l2v_sys_topwwwjetsdata_8TeV_{codename[0]}       lnN     -               """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)

    for i in range(0,Nbkg_int):
        if (i==1): card += """\t\t\t\t1.2""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
        else: card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)

    for j in range(0,Nbkg_int):
        for i in range(0,len(background_shapeSyst[j])):
            #change name here:
#            if (isItCorrelated(background_shapeSyst[j][i])):
#                name_forCorr=background_shapeSyst[j][i]
#                name_forCorr=name_forCorr.replace('zll_','')
#                name_forCorr=name_forCorr.replace('wz3lnu_','')
#                card += """
#{background_shapeSyst} shape  """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i,background_shapeSyst=name_forCorr)
            if not(isItCorrelated(background_shapeSyst[j][i])):
                card += """
{background_shapeSyst} shape1  """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i,background_shapeSyst=background_shapeSyst[j][i])
#                card += """
#{background_shapeSyst} shape  """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i,background_shapeSyst=background_shapeSyst[j][i])
                for k in range(0,j+1):
                    card += """-\t\t\t\t""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i,background_shapeSyst=background_shapeSyst[j][i])
                card += """1.0 """.format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i,background_shapeSyst=background_shapeSyst[j][i])
                for k in range(1,Nbkg_int-j):
                    card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i,background_shapeSyst=background_shapeSyst[j][i])

    card += """
QCDJeT_aTG           lnN     1.12          """
    for i in range(0,Nbkg_int):
        card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs,i=i)

    card += """
QCDscale_VV           lnN     -          """
    for i in range(0,Nbkg_int):
        if (i==2): card += """\t\t\t\t1.054""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
        else: card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)

    card += """
pdf_VV           lnN     1.058          """
    for i in range(0,Nbkg_int):
        if (i==2): card += """\t\t\t\t1.042""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
        else: card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
        
    if (doSignalShape_unc):
        for i in range(0,len(signal_shapeSyst)):
            #change name here:
            if (isItCorrelated(str(signal_shapeSyst[i]))):
                name_forCorr=str(signal_shapeSyst[i])
                name_forCorr=name_forCorr.replace('zz2l2nu_','')
                card += """
{signal_shapeSyst}        shape1  1.0          """.format(signal_shapeSyst=name_forCorr)
                for j in range(0,Nbkg_int):
                    isItCorr=False
                    for k in range(0,len(background_shapeSyst[j])):
                        if (name_forCorr in background_shapeSyst[j][k]):
                            isItCorr=True
                    if (isItCorr):
                        card += """\t\t\t\t1.0""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)
                    else:
                        card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)

            else:
                card += """
{signal_shapeSyst}        shape1  1.0          """.format(signal_shapeSyst=signal_shapeSyst[i])    
                for i in range(0,Nbkg_int):
                    card += """\t\t\t\t-""".format(codename=codename,norm_sig_sm=norm_sig_sm,norm_bkg=norm_bkg,norm_obs=norm_obs)

    print card

    cardfile = open('wv_semil_%sboosted.txt'%(codename),'w')
    cardfile.write(card)
    cardfile.close

