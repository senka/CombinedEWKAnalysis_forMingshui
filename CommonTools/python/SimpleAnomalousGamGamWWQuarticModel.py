from CombinedEWKAnalysis.CommonTools.AnomalousCouplingModel import *

class SimpleAnomalousGamGamWWQuarticModel(AnomalousCouplingModel):
    def __init__(self):
        AnomalousCouplingModel.__init__(self)
        self.processes = ['aaWWBSM']
        self.channels  = ['emu']
        self.pois    = ['a0W','aCW']        
        self.anomCoupSearchWindows = {'a0W':['-1e-3','1e-3'],
                                      'aCW':['-5e-3','5e-3']     }
        
        self.verbose = False

    def buildScaling(self,process,channel):
        scalerName = 'aaWWBSM'
        self.modelBuilder.factory_('expr::Scaling_%s("((7.85157e+07)'\
                                   '  * @0 * @0) + (2 * @0 * @1 *'\
                                   ' (1.94746e+07)) + ((5.90419e+06)'\
                                   ' * @1 * @1) +  0.0", a0W, aCW)'%scalerName)
        return scalerName
        

myModel = SimpleAnomalousGamGamWWQuarticModel()
