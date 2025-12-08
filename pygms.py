import json, bisect, os, sys, copy
sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))
import numpy as np
import pandas as pd
import USGSHazardGMM, SignificantDurationModel, CorrelationModel
# added
from util import *

LOCAL_IM_GMPE = {"DS575": ["Bommer, Stafford & Alarcon (2009)", "Afshari & Stewart (2016)"],
                 "DS595": ["Bommer, Stafford & Alarcon (2009)", "Afshari & Stewart (2016)"],
				 "DS2080": ["Afshari & Stewart (2016)"]}

AVAIL_IMG = ['SA','DS575','DS595']


class TargetIntensityMeasure:
    def __init__(self,tgt_im_config):
        """
        Input:
        - tgt_im_config: configuration dictionary with following attributes
        -- TargetType: target type, current options: 
           'CMS' (conditional mean spectrum), 'CS' (conditional spectrum), 'CSD' (conditional spectrum & duration)
        -- IntensityMeasureType: a list of intensity measure types, e.g., ['SA','DS575']
        -- Hazard: a dictionary including hazard information
        -- GroundMotionModel: a dictionary defining the ground motion model for individual intensity measures
        """
        # target type
        self.tgt_type = tgt_im_config.get('TargetType',None)
        if self.tgt_type is None:
            print('TargetIntensityMeasure.__init__: please define TargetType.')
            return
        elif self.tgt_type not in ['CMS','CS','CSD']:
            print('TargetIntensityMeasure.__init__: current available target type: CMS, CS, CSD.')
            return
        else:
            # configure conditional intensity measure type
            cond_im = tgt_im_config.get('ConditionalIntensityMeasure',None)
            if cond_im is None:
                print('TargetIntensityMeasure.__init__: please define ConditionalIntensityMeasure.')
                return
            else:
                cond_imt = list(cond_im.keys())
                if len(cond_imt) > 1:
                    print('TargetIntensityMeasure.__init__: only one conditional intensity measure currently supported.')
                    return
                else:
                    self.cond_imt = cond_imt[0]
                if self.cond_imt not in AVAIL_IMG:
                    print('TargetIntensityMeasure.__init__: current supported IM: {}.'.format(AVAIL_IMG))
                    return
                # conditional period
                self.cond_period = cond_im.get(self.cond_imt).get('Period',None)
                # conditional value
                self.cond_imv = cond_im.get(self.cond_imt).get('Value')

        # kz-250722: csv file input for periods
        self.filedirs = tgt_im_config.get('Directories',None)
        if self.filedirs is not None:
            self.input_dir = self.filedirs.get('Input',None)
            self.output_dir = self.filedirs.get('Output',None)
        else:
            self.input_dir = None
            self.output_dir = None
        # check input directory
        if self.input_dir is not None:
            if os.path.isdir(self.input_dir):
                pass
            else:
                print('TargetIntensityMeasure.__init__: ERROR - the TargetIntensityMeasure input directory {} not found.'.format(self.input_dir))
                return
        if self.output_dir is not None:
            if os.path.isdir(self.output_dir):
                os.makedirs(self.output_dir)
                print('TargetIntensityMeasure.__init__: MESSAGE - the TargetIntensityMeasure output directory {} created.'.format(self.output_dir))
                
        # configure intensity measure type
        self.imt = tgt_im_config.get('IntensityMeasureType',None)
        if self.imt is None:
            print('TargetIntensityMeasure.__init__: please define IntensityMeasureType.')
            return
        self.__imt_config()

        # configure hazard information
        self.hazard_info = tgt_im_config.get('Hazard',None)
        if self.hazard_info is None:
            print('TargetIntensityMeasure.__init__: please define Hazard.')
            return
        self.__hazard_config()

        # configure ground motion model:
        self.gmm_info = tgt_im_config.get('GroundMotionModel',None)
        if self.gmm_info is None:
            print('TargetIntensityMeasure.__init__: please define GroundMotionModel.')
            return
        self.__gmm_config()

        # configure intensity measure calculators
        self.set_im_calculator()
    

    def __imt_config(self):
        """
        Configure intensity measure types
        """
        # initialization
        self.ims = []
        # intensity measure groups
        self.im_groups = list(self.imt.keys())
        self.num_img = len(self.im_groups)
        if self.num_img < 1:
            print('TargetIntensityMeasure.__imt_config: no IMGroups found.')
            return
        # set individual intensity measures
        self.sa_periods = []
        for cur_img in self.im_groups:
            if cur_img not in AVAIL_IMG:
                print('TargetIntensityMeasure.__imt_config: {} is not in the current supported IMs: {}.'.format(cur_img, AVAIL_IMG))
                return
            if cur_img == 'SA':
                #sa_periods = self.imt.get(cur_img).get('Periods')
                container_periods = self.imt.get(cur_img).get('Periods')
                # kz-250722: csv input option for periods
                if type(container_periods) is list:
                    # directly defining the list of periods
                    self.sa_periods = container_periods
                elif type(container_periods) is str:
                    # pointing to a text file of periods
                    period_filepath = os.path.join(self.input_dir,container_periods)
                    if os.path.isfile(period_filepath):
                        # load the periods
                        self.sa_periods = load_user_defined_periods(period_filepath)
                    else:
                        print('TargetIntensityMeasure.__imt_config: ERROR - the user-defined period file {} not found.'.format(period_filepath))
                        return
                # add the conditional period if not already in the sa_periods
                if (self.cond_period) and (self.cond_period not in self.sa_periods):
                    bisect.insort(self.sa_periods,self.cond_period)
                if len(self.sa_periods) < 1:
                    print('TargetIntensityMeasure.__imt_config: please define Sa Periods.')
                    return
                self.ims = self.ims+['SA({})'.format(str(x).rstrip('0').rstrip('.')) for x in self.sa_periods]
            elif cur_img.startswith('DS'):
                self.ims.append(cur_img+'H')
            else:
                pass
    

    def __hazard_config(self):
        """
        Configure hazard information
        """
        # earthquake magnitude
        self.magnitude = self.hazard_info.get('Magnitude',None)
        # site to rupture distances
        self.rrup = self.hazard_info.get('Rrup',None)
        self.rjb = self.hazard_info.get('Rjb',None)
        self.rx = self.hazard_info.get('Rx',None)
        # vs30
        self.vs30 = self.hazard_info.get('Vs30',None)
        # check
        if not all([self.magnitude, self.vs30]) or not any([self.rrup, self.rjb, self.rx]):
            print('TargetIntensityMeasure: please provide Magintude, Vs30 and one distance at least.')
            return
        # other parameters
        self.dip = self.hazard_info.get('Dip',None)
        self.rake = self.hazard_info.get('Rake',None)
        self.width = self.hazard_info.get('Width',None)
        self.ztop = self.hazard_info.get('Ztop',None)
        self.zhyp = self.hazard_info.get('Zhyp',None)
        self.z1p0 = self.hazard_info.get('Z1p0',None)
        self.z2p5 = self.hazard_info.get('Z2p5',None)
    

    def __gmm_config(self):
        """
        Configure ground motion model
        """
        # set up ground motion prediction equations
        self.gmm = dict()
        for cur_im, cur_prop in self.gmm_info.items():
            if cur_prop is None:
                print('TargetIntensityMeasure.__gmm_config: cannot find GMPE {} not supported for {} yet.'.format(cur_prop,cur_im))
                return
            self.gmm.update({cur_im: cur_prop})
    

    def set_im_calculator(self):
        # loop over all imt
        self.im_calculators = dict()
        for cur_imt, cur_prop in self.imt.items():
            #cur_periods = cur_prop.get('Periods', None)
            if cur_imt == "SA":
                cur_periods = self.sa_periods
            else:
                cur_periods = cur_prop.get('Periods', None)
            cur_gmm = self.gmm.get(cur_imt, None)
            if cur_gmm is None:
                print('SiteData.compute_im: cannot find the ground motion model for {}'.format(cur_imt))
                return
            if cur_gmm == "Bommer, Stafford & Alarcon (2009)":
                cur_calculator = SignificantDurationModel.bommer_stafford_alarcon_ds_2009
            elif cur_gmm == "Afshari & Stewart (2016)":
                cur_calculator = SignificantDurationModel.afshari_stewart_ds_2016
            else:
                cur_calculator = USGSHazardGMM.USGS_GMM(cur_gmm, self.magnitude, self.vs30, imt=cur_imt, periods=cur_periods, rrup=self.rrup, 
                                                        rjb=self.rjb, rx=self.rx, dip=self.dip, width=self.width, ztop=self.ztop, 
                                                        zhyp=self.zhyp, z1p0=self.z1p0, z2p5=self.z2p5, rake=self.rake)
            # collect
            self.im_calculators.update({cur_imt: cur_calculator})
    

    def run_im_calculator(self,condition_flag=True):
        # initialize mean and standard deviation list
        im_mean = []
        im_std = []
        # loop over all intensity measure types
        for cur_imt, cur_prop in self.imt.items():
            cur_calculator = self.im_calculators.get(cur_imt)
            cur_gmm = self.gmm.get(cur_imt, None)
            if cur_gmm == "Bommer, Stafford & Alarcon (2009)":
                theta, beta_tot, beta_tau, beta_phi = cur_calculator(magnitude=self.magnitude, distance=self.rrup, 
                                                                     vs30=self.vs30, duration_type=cur_imt+'H')
                # collect
                im_mean.append(float(theta))
                im_std.append(float(beta_tot))
            elif cur_gmm == "Afshari & Stewart (2016)":
                theta, beta_tot, beta_tau, beta_phi = cur_calculator(magnitude=self.magnitude, distance=self.rrup, 
                                                                     vs30=self.vs30, duration_type=cur_imt+'H')
                # collect
                im_mean.append(float(theta))
                im_std.append(float(beta_tot))
            else:
                cur_calculator.fetch_data()
                theta = np.log(cur_calculator.im['Median']).tolist()
                beta_tot = cur_calculator.im['Dispersion']
                # collect
                im_mean = im_mean+theta
                im_std = im_std+beta_tot
        # correlation coefficient
        im_corr = np.zeros((len(self.ims),len(self.ims)))
        for j, im1 in enumerate(self.ims):
            for k, im2 in enumerate(self.ims):
                if j > k:
                    continue
                im_corr[j,k] = CorrelationModel.baker_bradley_correlation_2017(im1=im1, im2=im2)
                im_corr[k,j] = im_corr[j,k]
        # evaluate the conditional intensity measures
        if condition_flag:
            if self.cond_imt.startswith('DS'):
                self.cond_im_idx = self.ims.index(self.cond_imt+'H')
            else:
                self.cond_im_idx = self.ims.index(self.cond_imt+'({})'.format(str(self.cond_period).rstrip('0').rstrip('.')))
            cmim, cov_cond = conditional_intensity_measure(im_mean,im_std,im_corr.tolist(),self.cond_im_idx,self.cond_imv)
        else:
            cmim = []
            cov_cond = []

        # collect
        self.im_target = {
            'Median': np.exp(im_mean).tolist(),
            'StandardDev': im_std, 
            'Correlation': im_corr.tolist(),
            'ConditionalMean': np.exp(cmim).tolist(), 
            'ConditionalCov': cov_cond.tolist()
        }


class GroundMotionSelection:
    def __init__(self,gms_config):
        """
        Input:
        - gms_config: configuration dictionary for ground motion selection:
        -- RecordNum: number of ground motion records
        -- Scaling: min. and max. scaling factors
        -- Database: ground motion database directory
        -- ErrorWeight: error weight of mean, standard deviation, and different intensity measures 
        """
        # number of records
        self.num_rec = gms_config.get('RecordNum')
        # scaling factors (default 0.01~100)
        self.scaling_factor_min = gms_config.get('Scaling').get('Min',0.01)
        self.scaling_factor_max = gms_config.get('Scaling').get('Max',100.0)
        # added
        self.magnitude = gms_config.get('Hazard').get('Magnitude')
        self.vs30 = gms_config.get('Hazard').get('Vs30')
        self.distance = gms_config.get('Hazard').get('Rrup')
        #
        self.magnitude_min = gms_config.get('Filtering').get('Magnitude_min',0.0)
        self.magnitude_max = gms_config.get('Filtering').get('Magnitude_max',10)
        self.vs30_min = gms_config.get('Filtering').get('Vs30_min',0.0)
        self.vs30_max = gms_config.get('Filtering').get('Vs30_max',10000)
        self.distance_min = gms_config.get('Filtering').get('Rrup_min',0)
        self.distance_max = gms_config.get('Filtering').get('Rrup_max',10000)
        # kz-250722: adding a flag to control histogram plot (default is true)
        self.filtering_hist = gms_config.get('Filtering').get('PlotHistogram',True)
        # kz-250722: moving the random seed inside the selection class (default is None)
        self.random_seed = gms_config.get('RandomSeed',None)
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
        # end of addition
        # ground motion database directiory
        self.gmdb_file = gms_config.get('Database',None)
        if self.gmdb_file is None:
            self.gmdb_file = {
                "SA": os.path.join(os.path.dirname(__file__),'gmdb/nga_west_psa.csv'),
                "Periods": os.path.join(os.path.dirname(__file__),'gmdb/nga_west_period.csv'),
                "DS575": os.path.join(os.path.dirname(__file__),'gmdb/nga_west_ds575.csv'),
                "DS595": os.path.join(os.path.dirname(__file__),'gmdb/nga_west_ds595.csv'),
                "Filename": os.path.join(os.path.dirname(__file__),'gmdb/nga_west_filename.csv'),
                # added
                "Magnitude": os.path.join(os.path.dirname(__file__),'gmdb/nga_west_Earthquake Magnitude.csv'),
                "Vs30": os.path.join(os.path.dirname(__file__),'gmdb/nga_west_Preferred Vs30 (m-s).csv'),
                "Distance": os.path.join(os.path.dirname(__file__),'gmdb/nga_west_EpiD (km).csv'),
                # end of addition
            }
        elif self.gmdb_file == 'NGA-West2':
            self.gmdb_file = {
                "SA": [os.path.join(os.path.dirname(__file__),'gmdb/nga_west2/psa_p1.csv'),
                       os.path.join(os.path.dirname(__file__),'gmdb/nga_west2/psa_p2.csv'),
                       os.path.join(os.path.dirname(__file__),'gmdb/nga_west2/psa_p3.csv')],
                "Periods": os.path.join(os.path.dirname(__file__),'gmdb/nga_west2/period.csv'),
                "DS575": os.path.join(os.path.dirname(__file__),'gmdb/nga_west2/ds575.csv'),
                "DS595": os.path.join(os.path.dirname(__file__),'gmdb/nga_west2/ds595.csv'),
                "Filename": os.path.join(os.path.dirname(__file__),'gmdb/nga_west2/rsn.csv'),
                # added
                "Magnitude": os.path.join(os.path.dirname(__file__),'gmdb/nga_west2/Earthquake Magnitude.csv'),
                "Vs30": os.path.join(os.path.dirname(__file__),'gmdb/nga_west2/Preferred Vs30 (m-s).csv'),
                "Distance": os.path.join(os.path.dirname(__file__),'gmdb/nga_west2/EpiD (km).csv')
            }
        self.gmdb = dict()
        self.gmdb_status = False
        # error weight
        self.err_weight = gms_config.get('ErrorWeight',None)
        if self.err_weight is None:
            self.err_weight = {
                "SA": [0.5, 0.5],
                "DS575": [0.5, 0.5]
            }
        self.tgt_type = gms_config.get('TargetType')
    

    def __load_gmdb(self, ims, cond_imv, cond_idx):
        # load the ground motion database
        for cur_key, cur_db_filename in self.gmdb_file.items():
            if isinstance(cur_db_filename,list):
                cur_db = pd.DataFrame()
                for cur_file in cur_db_filename:
                    tmp = pd.read_csv(cur_file,index_col=False,header=None)
                    cur_db = pd.concat([cur_db, tmp], ignore_index=True)
            else:
                cur_db = pd.read_csv(cur_db_filename,index_col=False,header=None)
            self.gmdb.update({cur_key: cur_db})
        self.gmdb_periods = self.gmdb.get('Periods').to_numpy().flatten().tolist()
        self.gmdb_imv = self.gmdb['SA'].to_numpy()

        # added
        self.gmdb_magnitude = self.gmdb['Magnitude'].to_numpy()
        self.gmdb_vs30 = self.gmdb['Vs30'].to_numpy()
        self.gmdb_distance = self.gmdb['Distance'].to_numpy()
        # end of addition

        # scaling factor check
        cond_ims = ims[cond_idx]
        if cond_ims.startswith('SA'):
            gmdb_cond_idx = self.gmdb_periods.index(float(cond_ims[3:-1]))
        self.gm_scaling = [float(cond_imv/x) for x in self.gmdb_imv[:,gmdb_cond_idx]]
        self.gm_flag = []
        
        #added (histograms)
        if self.filtering_hist:
            plot_filter_effect(self.gm_scaling, self.scaling_factor_min, self.scaling_factor_max, "Scaling Factor", None)
            plot_filter_effect(self.gmdb_magnitude, self.magnitude_min, self.magnitude_max, "Magnitude", self.magnitude)
            plot_filter_effect(self.gmdb_vs30, self.vs30_min, self.vs30_max, "Vs30", self.vs30)
            plot_filter_effect(self.gmdb_distance, self.distance_min, self.distance_max, "Rrup (Distance)", self.distance)
        # end of addition

        # modified
        # modification to filter with magnitude, Vs30 and Rrup, and to monitor filter effect
        print(f"Total initial samples: {len(self.gm_scaling)}")
        count_after_scaling = 0
        count_after_magnitude = 0
        count_after_vs30 = 0
        count_after_distance = 0
        for gmid,cur_scaling in enumerate(self.gm_scaling):
            cur_magnitude=self.gmdb_magnitude[gmid]
            cur_vs30=self.gmdb_vs30[gmid]
            cur_distance=self.gmdb_distance[gmid]
            if self.scaling_factor_min <= cur_scaling <= self.scaling_factor_max:
                count_after_scaling += 1
                if self.magnitude_min <= cur_magnitude <= self.magnitude_max:
                    count_after_magnitude += 1
                    if self.vs30_min <= cur_vs30 <= self.vs30_max:
                        count_after_vs30 += 1
                        if self.distance_min <= cur_distance <= self.distance_max:
                            count_after_distance += 1
                            self.gm_flag.append(gmid)
        self.gm_sf = np.array([self.gm_scaling[x] for x in self.gm_flag])
        self.gmdb_imv = np.log(self.gmdb_imv[self.gm_flag,:]*self.gm_sf[:,np.newaxis])
        # end of modification
        # added
        print("Filtering summary:")
        print(f"After scaling filter: {count_after_scaling}")
        print(f"After magnitude filter: {count_after_magnitude}")
        print(f"After Vs30 filter: {count_after_vs30}")
        print(f"After distance filter: {count_after_distance}")
        print(f"Final accepted samples: {len(self.gmdb_imv)}")
        # end of addition
        for cur_im in ims:
            if cur_im.startswith('DS575'):
                tmp = np.log(self.gmdb['DS575'].to_numpy())
                self.gmdb_imv = np.concatenate((self.gmdb_imv, tmp[self.gm_flag,:]), axis=1)
            elif cur_im.startswith('DS595'):
                tmp = np.log(self.gmdb['DS595'].to_numpy())
                self.gmdb_imv = np.concatenate((self.gmdb_imv, tmp[self.gm_flag,:]), axis=1)
            else:
                pass    
        self.gmdb_status = True

    
    def __unload_gmdb(self):
        # unload the ground motion database
        self.gmdb = dict()
        self.gmdb_status = False

    
    def generate_pseudo_im(self,dict_tgt_im):
        # generate pseudo intensity measures
        self.tgt_mean = np.log(dict_tgt_im.get('ConditionalMean'))
        self.tgt_cov = np.array(dict_tgt_im.get('ConditionalCov'))
        self.pseudo_im_rlz = np.random.multivariate_normal(self.tgt_mean, self.tgt_cov, self.num_rec)


    def initial_scanning(self, ims, cond_imv, cond_idx):
        """
        Search the best match for each realization of pseudo intensity measures
        """
        # load the ground motion database if not yet
        if not self.gmdb_status:
            self.__load_gmdb(ims, cond_imv, cond_idx)
        self.selection_periods = [float(x[3:-1]) for x in ims if x.startswith('SA')]
        self.ims_idx = [self.gmdb_periods.index(x) for x in self.selection_periods]
        for cur_im in ims:
            if cur_im.startswith('DS575'):
                self.ims_idx.append(len(self.gmdb_periods))
            elif cur_im.startswith('DS595'):
                self.ims_idx.append(len(self.gmdb_periods.shape)+1)
        # selected_id
        self.selected_id = []
        for i in range(self.num_rec):
            cur_tgt = self.pseudo_im_rlz[i,:]
            ind_err = np.abs(self.gmdb_imv[:,self.ims_idx]-cur_tgt[np.newaxis,:])
            if 'DS575H' in ims:
                cur_err = np.mean(ind_err[:,:-1],axis=1)*self.err_weight['SA'][0]+ind_err[:,-1]*self.err_weight['DS575'][0]
            elif 'DS595H' in ims:
                cur_err = np.mean(ind_err[:,:-1],axis=1)*self.err_weight['SA'][0]+ind_err[:,-1]*self.err_weight['DS595'][0]
            else:
                cur_err = np.mean(ind_err[:,:-1],axis=1)
            cur_sort = np.argsort(cur_err)
            for cur_id  in cur_sort:
                if cur_id not in self.selected_id:
                    self.selected_id.append(int(cur_id))
                    break
                else:
                    pass

    def optimize_selection(self, ims, num_greedy_rounds=2):
        """
        Greedy algorithm search to minimize the combined error
        (fast version, unsafe updates, optimized for performance)

        Changes vs optimize_selection2:
        - Precomputes gmdb_subset and target std once (avoids repeated slicing and diag calls)
        - Precomputes weights once before loops (avoids repeated if checks)
        - More in-place trial updates for speed (less copying)
        """
        cur_selected_id = copy.deepcopy(self.selected_id)
        # Precompute subset matrix once for faster access (new in v3)
        gmdb_subset = self.gmdb_imv[:, self.ims_idx]  # shape (n_motions, n_ims), avoids long re-slicing
        tgt_std = np.sqrt(np.diag(self.tgt_cov))      # precompute target std once (new in v3)
        # Precompute weights once before loops (new in v3)
        if self.tgt_type == 'CSD':
            if 'DS575H' in ims:
                w_sa_mean, w_ds_mean = self.err_weight['SA'][0], self.err_weight['DS575'][0]
                w_sa_std,  w_ds_std  = self.err_weight['SA'][1], self.err_weight['DS575'][1]
            elif 'DS595H' in ims:
                w_sa_mean, w_ds_mean = self.err_weight['SA'][0], self.err_weight['DS595'][0]
                w_sa_std,  w_ds_std  = self.err_weight['SA'][1], self.err_weight['DS595'][1]
            else:
                w_sa_mean = self.err_weight['SA'][0]
                w_sa_std = self.err_weight['SA'][1]
                w_ds_mean = w_ds_std = 0.0
        else:
            w_sa_mean = self.err_weight['SA'][0]
            w_sa_std = self.err_weight['SA'][1]
            w_ds_mean = w_ds_std = 0.0
        # Initial error calculation (same logic as v2)
        selected_matrix = gmdb_subset[cur_selected_id]
        mean_vec = np.mean(selected_matrix, axis=0)
        std_vec = np.std(selected_matrix, axis=0)
        err_mean = np.square(mean_vec - self.tgt_mean)
        err_std = np.square(std_vec - tgt_std)
        err_tot = (np.mean(err_mean[:-1]) * w_sa_mean + err_mean[-1] * w_ds_mean + np.mean(err_std[:-1]) * w_sa_std + err_std[-1] * w_ds_std)
        print(f'GroundMotionSelection.optimize_selection: current total error: {err_tot:.6f}.')
        for i in range(num_greedy_rounds):
            for cur_gmseat in range(self.num_rec):
                print(f'GroundMotionSelection.optimize_selection: {100*(self.num_rec*i+(cur_gmseat+1))/self.num_rec/num_greedy_rounds:.2f}%.')
                old_idx = cur_selected_id[cur_gmseat]
                # Precompute sums and sumsq once per seat (as in v2)
                selected_matrix = gmdb_subset[cur_selected_id]
                sum_selected = np.sum(selected_matrix, axis=0)
                sumsq_selected = np.sum(np.square(selected_matrix), axis=0)
                old_row = gmdb_subset[old_idx]
                for new_gmidx in range(len(self.gm_flag)):
                    if new_gmidx in cur_selected_id:
                        continue  # skip if already selected
                    new_row = gmdb_subset[new_gmidx]
                    cur_selected_id[cur_gmseat] = new_gmidx  # trial replace (more in-place in v3)
                    # Update sums and means without full recompute (same as v2)
                    sum_adj = sum_selected - old_row + new_row
                    mean_adj = sum_adj / self.num_rec
                    err_mean = np.square(mean_adj - self.tgt_mean)
                    # Update std using adjusted sumsq (same as v2)
                    sumsq_adj = sumsq_selected - old_row ** 2 + new_row ** 2
                    std_sq = sumsq_adj / self.num_rec - mean_adj ** 2
                    std_adj = np.sqrt(np.maximum(std_sq, 0.0))
                    err_std = np.square(std_adj - tgt_std)
                    # Combined error (same formula as v2, but weights precomputed)
                    err_new = (np.mean(err_mean[:-1]) * w_sa_mean + err_mean[-1] * w_ds_mean + np.mean(err_std[:-1]) * w_sa_std + err_std[-1] * w_ds_std)
                    if err_new < err_tot:
                        self.selected_id = copy.deepcopy(cur_selected_id)
                        err_tot = err_new
                # MD Correction
                cur_selected_id = copy.deepcopy(self.selected_id)
                print(f'GroundMotionSelection.optimize_selection: current total error: {err_tot:.6f}.')


    def export_selection(self, ims, output_dir=None):
        if output_dir is None:
            output_dir = os.path.dirname(__file__)
            output_path = os.path.join(output_dir,'GroundMotionSelectionResults.csv')
        # ground motion names
        self.selected_gm_names = self.gmdb['Filename'].iloc[[self.gm_flag[x] for x in self.selected_id],0].to_list()
        # intensity measures
        self.selected_im = self.gmdb['SA'].iloc[[self.gm_flag[x] for x in self.selected_id],:].to_numpy()
        if self.tgt_type in ['CSD']:
            if 'DS575H' in ims:
                self.selected_im = np.concatenate((self.selected_im,self.gmdb['DS575'].iloc[[self.gm_flag[x] for x in self.selected_id],:].to_numpy()),axis=1)
            elif 'DS595H' in ims:
                self.selected_im = np.concatenate((self.selected_im,self.gmdb['DS595'].iloc[[self.gm_flag[x] for x in self.selected_id],:].to_numpy()),axis=1)
        # scaling factor
        self.selected_sf = [self.gm_scaling[x] for x in [self.gm_flag[x] for x in self.selected_id]]
        # output df
        dict_output = {
            'GroundMotion': self.selected_gm_names
        }
        for i,cur_period in enumerate(self.gmdb_periods):
            dict_output.update({
                'SA({})'.format(cur_period): self.selected_im[:,i].tolist()
            })
        if self.tgt_type in ['CSD']:
            if 'DS575H' in ims:
                dict_output.update({
                    'DS575': self.selected_im[:,-1].tolist()
                })
            elif 'DS595H' in ims:
                dict_output.update({
                    'DS595': self.selected_im[:,-1].tolist()
                })
        dict_output.update({
            'ScalingFactor': self.selected_sf
        })
        df_output = pd.DataFrame.from_dict(dict_output)
        # save
        df_output.to_csv(output_path)


def conditional_intensity_measure(im_mean,im_std,im_corr,im_cond_idx,im_cond):
    """
    Input:
    - im_mean: the mean value of logarithmic intensity measures
    - im_std: the total standard deviation of logarithmic intensity measures
    - im_corr: the correlation coefficient matrix of logarithmic intensity measures
    - im_cond_idx: the index of conditional intensity measure (e.g., Sa(1.0s) is the 20th number in the im_mean)
    - im_cond: the value of the conditional intensity measure (e.g., would like to target the Sa(1.0s) to 0.9g)
    Output:
    - cmim: conditional mean intensity measures (in log scale)
    - cov_cond: covariance matrix conditional on the given intensity measure
    """
    # conditional mean spectrum
    eps_bar = (np.log(im_cond)-im_mean[im_cond_idx])/im_std[im_cond_idx]
    cmim = np.array(im_mean)+np.multiply(np.array(im_std),np.array(im_corr[im_cond_idx]))*eps_bar
    # covariance matrix
    cov_matrix = np.diag(im_std).dot(np.array(im_corr).dot(np.diag(im_std)))
    # partition matrix to submatrices
    cov_11 = cov_matrix[im_cond_idx:im_cond_idx+1,im_cond_idx:im_cond_idx+1]
    cov_12 = np.delete(cov_matrix,im_cond_idx,axis=1)[im_cond_idx:im_cond_idx+1,:]
    cov_21 = np.delete(cov_matrix,im_cond_idx,axis=0)[:,im_cond_idx:im_cond_idx+1]
    cov_22 = np.delete(np.delete(cov_matrix,im_cond_idx,axis=0),im_cond_idx,axis=1)
    # conditional matrix
    cov_cond = cov_22 - np.dot(cov_21, np.linalg.inv(cov_11)).dot(cov_12)
    # insert the conditional row and column
    cov_cond = np.insert(cov_cond,im_cond_idx,0,axis=0)
    cov_cond = np.insert(cov_cond,im_cond_idx,0,axis=1)
    if not np.allclose(cov_cond, cov_cond.T):
        print('conditional_intensity_measure: covariance matrix is not symmetric...')
    eigenvalues = np.linalg.eigvalsh(cov_cond)
    if np.all(eigenvalues > -1e-8):
        pass
    else:
        print('conditional_intensity_measure: covariance matrix is not positive semidefinite...')
    # return
    return cmim, cov_cond

# adding back the run_selection in pygms for batch run (kz)
def run_selection(gms_config):
    with open(gms_config,'r') as f:
        job_config = json.load(f)

    tgt_config = job_config.get('TargetIntensityMeasure')
    my_tgt = TargetIntensityMeasure(tgt_config)
    my_tgt.run_im_calculator()

    gms_config = job_config.get('GroundMotionSelection')
    gms_config.update({'TargetType':tgt_config.get('TargetType')})
    gms_config['Hazard']=tgt_config.get('Hazard')

    my_gms = GroundMotionSelection(gms_config)
    my_gms.generate_pseudo_im(my_tgt.im_target)
    my_gms.initial_scanning(my_tgt.ims,my_tgt.cond_imv,my_tgt.cond_im_idx)
    my_gms.optimize_selection(my_tgt.ims)

    my_gms.export_selection(my_tgt.ims)