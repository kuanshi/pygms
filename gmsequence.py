from pygms import *
from scipy.interpolate import RegularGridInterpolator
sys.path.append(os.path.join(os.path.dirname(__file__), "lib/sequence_corr"))

# default code for periods: -3: Ds595, -2: Ds575, -1: PGV, 0: PGA, >0: Sa(T)
PERIODS = np.array([-3,-2,-1,0]+[0.01+0.01*x for x in range(1000)])
# Available ground motion pair correlation models
with open('available_models.json','r') as f:
    dict_seq_corr_models = json.load(f)
LOCAL_CORR_MODEL = dict_seq_corr_models.get('CorrelationModels')
LOCAL_CORR_FILE = dict_seq_corr_models.get('ModelFiles')

class TargetIntensityMeasureSequence:
    """
    GroundMotionSequence: a sequence of ground motions at a target site location in an earthquake sequence
    """
    def __init__(self,gmseq_config):
        # number of events in the sequence (default is 1 - a single earthquake event)
        self.num_events = gmseq_config.get('NumberEvents',1)
        # target intensity measure sequence
        self.tgt_im_sequence = gmseq_config.get('TargetIntensityMeasureSequence',[])
        if len(self.tgt_im_sequence) != self.num_events:
            print('GroundMotionSequence.__init__: please check the target TargetIntensityMeasureSequence list whose number does not match NumberEvents.')
            return
        # set hazards
        self.__configure_hazards()
        # sequential-gm correlation model
        self.seq_corr = gmseq_config.get('SequenceCorrelation',dict())
        # set correlation models
        self.__configure_gmseq_correlation()
        # set target intensity measures
        self.__configure_tgt_intensity_measures()


    def __configure_hazards(self):
        """
        Configure hazard information
        """
        # earthquake magnitude
        self.magnitudes = [eq.get('Hazard').get('Magntidue') for eq in self.tgt_im_sequence]
        # site to rupture distances
        self.rrups = [eq.get('Hazard').get('Rrup',None) for eq in self.tgt_im_sequence]
        self.rjbs = [eq.get('Hazard').get('Rjb',None) for eq in self.tgt_im_sequence]
        self.rxs = [eq.get('Hazard').get('Rx',None) for eq in self.tgt_im_sequence]
        # vs30
        self.vs30s = [eq.get('Hazard').get('Vs30',None) for eq in self.tgt_im_sequence]
        # other parameters
        self.dips = [eq.get('Hazard').get('Dip',None) for eq in self.tgt_im_sequence]
        self.rakes = [eq.get('Hazard').get('Rake',None) for eq in self.tgt_im_sequence]
        self.widths = [eq.get('Hazard').get('Width',None) for eq in self.tgt_im_sequence]
        self.ztops = [eq.get('Hazard').get('Ztop',None) for eq in self.tgt_im_sequence]
        self.zhyps = [eq.get('Hazard').get('Zhyp',None) for eq in self.tgt_im_sequence]
        self.z1p0s = [eq.get('Hazard').get('Z1p0',None) for eq in self.tgt_im_sequence]
        self.z2p5s = [eq.get('Hazard').get('Z2p5',None) for eq in self.tgt_im_sequence]
        # CRjb matrix (in km)
        self.CRJBs = [eq.get('Hazard').get('CRJB',None) for eq in self.tgt_im_sequence]
        # time interval (in days)
        self.dTs = [eq.get('Hazard').get('dT',None)for eq in self.tgt_im_sequence]
        # azmuth angle
        self.dAs = [eq.get('Hazard').get('AzmuthAngle',None) for eq in self.tgt_im_sequence]


    def __configure_gmseq_correlation(self):
        # MS-AS correlation model type:
        self.corr_type_c1c2 = self.seq_corr.get('Class1-Class2',None)
        # between aftershocks:
        self.corr_type_c2c2 = self.seq_corr.get('Class2-Class2',None)
        self.corr_type_c21c21 = self.seq_corr.get('Class2.1-Class2.1',None)
        self.corr_type_c21c22 = self.seq_corr.get('Class2.1-Class2.2',None)


    def __configure_tgt_intensity_measures(self):
        # initialize individual event's target intensity measures
        self.indiv_intensity_measures = []
        for idx in range(self.num_events):
            self.indiv_intensity_measures.append(TargetIntensityMeasure(self.tgt_im_sequence[idx]))


    def compute_c1c2_corr(self,ctype=None):
        """
        Compute the full correlation coefficient matrix for class1-class2 ground motion pairs
        """
        self.corr_c1c2_periods = PERIODS
        if ctype is None:
            # uncorrelated assumption            
            self.corr_c1c2 = np.zeros((len(PERIODS),len(PERIODS)))
        elif ctype not in LOCAL_CORR_MODEL.get('Class1-Class2'):
            print('TargetIntensityMeasureSequence.compute_c1c2_corr: {} is not supported - please select one from {}.'.format(ctype,LOCAL_CORR_MODEL.get('Class1-Class2')))
            return
        else:
            # initialize
            self.corr_c1c2 = np.zeros((len(PERIODS),len(PERIODS)))
            # get the model files(s)
            cm_file = LOCAL_CORR_FILE.get(ctype)
            if cm_file.endswith('.txt'):
                # directly loading and get the discrete model first
                tmp = np.loadtxt(cm_file,delimiter=' ')
                tmp_periods = tmp[:,0]
                Tmin = np.min(tmp_periods)
                Tmax = np.max(tmp_periods)
                tmp_corr = tmp[:,1:]
                # find periods within the tmp_periods
                period_idx = [idx for idx,Ti in enumerate(self.corr_c1c2_periods) if Ti>=Tmin and Ti<=Tmax]
                # 2d-interpolate (in logT scale)
                cur_interp = RegularGridInterpolator((tmp_periods, tmp_periods), tmp_corr)
                Xi, Yi = np.meshgrid(self.corr_c1c2_periods[period_idx], self.corr_c1c2_periods[period_idx], indexing='ij')
                cur_locs = np.column_stack((Xi.ravel(),Yi.ravel()))
                self.corr_c1c2[np.ix_(period_idx,period_idx)] = cur_interp(cur_locs).reshape(len(period_idx),len(period_idx))
            elif cm_file.endswith('.py'):
                # a python function with input of (Ti,Tj)
                pass