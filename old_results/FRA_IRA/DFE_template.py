import dadi
import dadi.DFE as DFE
import pickle
import nlopt
from pathlib import Path
import matplotlib.pyplot as plt
from dadi import Numerics, Integration, PhiManip, Spectrum

# # Use dadi functions to parse the VCF file to generate a data dictionary.
# datafile = '/groups/rgutenk/emfonseca/projects/FishProject/Fish_DFE/Pbimac/Data/Pbimac.biallelic.missing75.nofixed.DP6.minQ30.no0or1homozy.recode.CDS.nonsynonymous.vcf.gz'
# dd = dadi.Misc.make_data_dict_vcf(datafile, '/groups/rgutenk/emfonseca/projects/FishProject/Fish_DFE/Pbimac/Data/Pbimac_poplist.txt')

# #Extract the SFS for the Yoruba (YRI) and Central European (CEU) population from that dictionary, with both projected down to 20. We project down for this example to make it doable on a laptop. For a real analysis we would probably not project so severely (ex. the 1000 Genomes Yoruba population's sample size can be as high as 216).
# pop_ids, ns = ['PbimacNS', 'PbimacS'], [20, 19]
# data_fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized=False)
# data_fs.mask[0,1] = True
# data_fs.mask[1,0] = True

# # Save our extracted spectrum to disk.
# data_fs.to_file('/groups/rgutenk/emfonseca/projects/FishProject/Fish_DFE/Pbimac/Data/Pbimac.biallelic.missing75.nofixed.DP6.minQ30.no0or1homozy.recode.CDS.nonsynonymous.fs')

# #We can see how much data is in the SFS. As a rule of thumb we want to maximize this number. We can do this by trying different sample sizes for our population(s).
# print(data_fs.S())


def main():
	# Load synonymous frequency spectrum
	data_fs = dadi.Spectrum.from_file('IRA_FRA_syn_unfolded.fs')

	# Retrive the sample sizes from the data
	ns = data_fs.sample_sizes

	# Define the synonymous theta you got from the demography fit
	theta0 = 26681.214011656422

	# Calculate the theta for the nonsynonymous data based on the ratio of nonsynonymous mutations to synonymous mutations. You will probably be get this ratio from the paper or calculate it, but typically it is larger grater than 1.
	theta_ns = theta0 * 2.4

	# Define the grid points based on the sample size. Because you are adding selection, you might want to increase the sizes of the grid points, as spectra with higher levels of selection are harder for dadi to calculate accurately
	pts_l = [max(ns)+140, max(ns)+150, max(ns)+160]

	# Get the selection version of the demographic model Depending if you want the version with selection being independent or shared between populations the model name is slightly different. Extrapolation happens in the process of making the cache, so you do not need to to wrap the extrapolation function. split_mig_sel is the version of the split_mig demographic model with independant selection

	demo_sel_model = DFE.DemogSelModels.IM_pre

	# Define the optimial parameters from the demography fits. If you used misidentification model, you can remove the misidentification parameter.
	demo_popt = [2.3892200757761812, 0.09865800019618395, 0.059546476404449086, 0.29863880148384436, 2.9956042751467993, 0.1765401593467742, 2.1701288461162407, 0.19321473756551447]
	# You could also loop through your results file and extract the parameters that way.

	# Generate cache. The gamma_bounds argument defines the range of the gamma distribution. The gamma_pts argument can be used to specify the number of selection coefficients that will be selected in that range to generate your cache. It is recommended to use gamma_bounds=[1e-4, 2000], gamma_pts=50 for either 1D or 2D cache generation on the HPC. If you want to generate the 2D cache (independent selection coefficients), use:
	# DFE.Cache2D
	# NOTE: When testing locally, having mp = True might cause a bug, so set it to mp=False until you start working on an HPC
	# cache2d = DFE.Cache2D(popt, ns, demo_sel_model, pts=pts_l, gamma_bounds=[1e-4, 2000], gamma_pts=50)

	cache2d = pickle.load(open('mmd_FRA_mmd_IRA_2d_cache.bpkl', 'rb'))

	# We can check if the cached spectra have any large negative values:
	if (cache2d.spectra<0).sum() > 0:
		print(
			'!!!WARNING!!!\nPotentially large negative values!\nMost negative value is: '+str(cache2d.spectra.min())+
			'\nIf negative values are very negative (<-0.001), rerun with larger values for pts_l'
			)

	# Define the DFE function used cache.integrate or DFE.mixture. Here we will use the cache.integration function
	dfe_func = cache2d.integrate

	# Define the selection distribution you want to use. If using the mixture model you will want to define the selection distribution for both 1D and 2D distributions.
	sele_dist2d = DFE.PDFs.biv_lognormal

	# Optimization for the DFE requires extra arguments that are not included in the optimizer function, so we need to define them ourselves.
	func_args = [sele_dist2d, theta_ns]

	# For the DFE.mixture function the argument would be along the lines of:
	# func_args = [cache1d, cache2d, sele_dist1d, sele_dist2d, theta_ns]

	# Choose starting parameters for inference. This is an example for the gamma distribution. Most importantly are the first two parameters shape and scale (also called alpha and beta).
	params = [0.2,0.2,0.9]

	# Define boundaries of optimization. It is a good idea to have boundaries for the DFE as the optimizer can take parameters to values that cause errors with calulating the spectrum. Ex. rho = 1 for independent population selection. If optimization runs up against boundaries, you can increase them
	lower_bounds = [-10, 0.01, 1e-1]
	upper_bounds = [10, 10, 0.99]

	for x in range(1):
		
		# For running on the HPC, it is a good idea to check if file is made before making it so that you don't overwrite other results
		# if x == 0:
			
		# 	print('Beginning optimization ************************************************')
			
		# 	fid = Path('/groups/rgutenk/emfonseca/projects/FishProject/Fish_DFE/Pbimac/Results/Model_fit/Pbimac_dfe_fits_sec_cont.txt')
		# 	fid = open(fid, "w")
		# 	silent = fid.write('log_likelihood\tmean\tvariance\tcorrelation_coefficient\ttheta\n')
		# 	fid.close()
		
		# Perturb parameters
		# Optimizers dadi uses are mostly deterministic
		# so you will want to randomize parameters for each optimization.
		# It is recommended to optimize at least 20 time if you only have
		# your local machin, 100 times if you have access to an HPC.
		# If you want a single script to do multiple runs, you will want to
		# start a for loop here
		p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,
									lower_bound=lower_bounds)
		
		# Run optimization
		# At the end of the optimization you will get the
		# optimal parameters and log-likelihood.
		# You can modify verbose to watch how the optimizer behaves,
		# what number you pass it how many evaluations are done
		# before the evaluation is printed.
		# For the DFE, because we calculated the theta_ns, we want to
		# set multinom=False.
		popt, ll_model = dadi.Inference.opt(p0, data_fs, dfe_func, pts=None, 
											func_args=func_args, 
											lower_bound=lower_bounds, 
											upper_bound=upper_bounds,
											maxeval=6000, multinom=False, verbose=0)
		
		# Generate DFE spectrum
		model_fs = dfe_func(popt, None, sele_dist2d, theta_ns, None)
		
		fid = open('test_biv_lognormal_dfe.txt','a')
		
		# Write results to fid
		res = [ll_model] + list(popt) + [theta_ns]
		silent = fid.write('\t'.join([str(ele) for ele in res])+'\n')
		fid.close()
		
		#print the iteration
		print('Finished optimization: {0}'.format(x+1))



	popt = [-0.37790407357018946, 0.18908242695311342, 0.9800403505231494]
	model_fs = dfe_func(popt, None, sele_dist2d, theta_ns, None)

	# Plot
	fig = plt.figure(219033)
	fig.clear()
	dadi.Plotting.plot_2d_comp_Poisson(model_fs, data_fs,vmin=1)
	# fig.savefig('')
if __name__ == '__main__':
	main()