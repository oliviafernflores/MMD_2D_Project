number="0051179"

cache1d_file="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/FRA_IRA/mmd_FRA_mmd_IRA_1d_cache.bpkl"
cache2d_file="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/FRA_IRA/mmd_FRA_mmd_IRA_2d_cache.bpkl"
boots_syn_dir="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_demography/IRA_FRA/bootstraps"
boots_nsyn_dir="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/GO_bootstraps/GO_${number}_bootstraps"
fs_file="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/GO_${number}_SFS/${number}_dict_Mmd_FRA_Mmd_IRA_SFS_nsyn_unfolded.fs"


# Symmetric Bivariate Lognormal
# Error: numpy.linalg.LinAlgError: Singular matrix
dfe_popt_file="mmd_FRA_mmd_IRA_2d_biv_lognormal_DFE_${number}.InferDFE.bestfits"
output_file="mmd_FRA_mmd_IRA_2d_biv_lognormal_DFE_${number}.InferDFE.godambe"
dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_lognormal --dfe-popt "$dfe_popt_file" --cache2d "$cache2d_file" --output "$output_file"


# Asymmetric Bivariate Lognormal
# Error: numpy.linalg.LinAlgError: Singular matrix
dfe_popt_file="mmd_FRA_mmd_IRA_2d_biv_lognormal_asymmetric_DFE_${number}.InferDFE.bestfits"
output_file="mmd_FRA_mmd_IRA_2d_biv_lognormal_asymmetric_DFE_${number}.InferDFE.godambe"
dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_lognormal --dfe-popt "$dfe_popt_file" --cache2d "$cache2d_file" --output "$output_file"


# Lognormal
# Error: None!
dfe_popt_file="mmd_FRA_mmd_IRA_1d_lognormal_DFE_${number}.InferDFE.bestfits"
output_file="mmd_FRA_mmd_IRA_1d_lognormal_DFE_${number}.InferDFE.godambe"
dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d lognormal --dfe-popt "$dfe_popt_file" --cache1d "$cache1d_file" --output "$output_file"


# Symmetric Bivariate Independent Gamma
# Error: None!
dfe_popt_file="mmd_FRA_mmd_IRA_2d_biv_ind_gamma_symmetric_DFE_${number}.InferDFE.bestfits"
output_file="mmd_FRA_mmd_IRA_2d_biv_ind_gamma_symmetric_DFE_${number}.InferDFE.godambe"
dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_ind_gamma --dfe-popt "$dfe_popt_file" --cache2d "$cache2d_file" --output "$output_file"


# Asymmetric Bivariate Independent Gamma
# Error: None!
dfe_popt_file="mmd_FRA_mmd_IRA_2d_biv_ind_gamma_asymmetric_DFE_${number}.InferDFE.bestfits"
output_file="mmd_FRA_mmd_IRA_2d_biv_ind_gamma_asymmetric_DFE_${number}.InferDFE.godambe"
dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_ind_gamma --dfe-popt "$dfe_popt_file" --cache2d "$cache2d_file" --output "$output_file"

# Gamma
# Error: None!
dfe_popt_file="mmd_FRA_mmd_IRA_2d_gamma_DFE_${number}.InferDFE.bestfits"
output_file="mmd_FRA_mmd_IRA_2d_gamma_DFE_${number}.InferDFE.godambe"
dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d gamma --dfe-popt "$dfe_popt_file" --cache1d "$cache1d_file" --output "$output_file"
