boots_nsyn_dir="GO_0009987_nsyn_bootstraps"
boots_syn_dir="IRA_FRA_syn_bootstraps"
fs_file="0009987_Mmd_FRA_Mmd_IRA_SFS_nsyn_unfolded.fs"
cache2d_file="mmd_FRA_mmd_IRA_2d_cache.bpkl"
cache1d_file="mmd_FRA_mmd_IRA_1d_cache.bpkl"

# Symmetric Bivariate Lognormal
# Error: numpy.linalg.LinAlgError: Singular matrix
dfe_popt_file="mmd_FRA_mmd_IRA_2d_biv_lognormal_DFE_0009987.InferDFE.bestfits"
output_file="mmd_FRA_mmd_IRA_2d_biv_lognormal_DFE_0009987.InferDFE.godambe"
dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_lognormal --dfe-popt "$dfe_popt_file" --cache2d "$cache2d_file" --output "$output_file"


# Asymmetric Bivariate Lognormal
# Error: numpy.linalg.LinAlgError: Singular matrix
dfe_popt_file="mmd_FRA_mmd_IRA_2d_biv_lognormal_asymmetric_DFE_0009987.InferDFE.bestfits"
output_file="mmd_FRA_mmd_IRA_2d_biv_lognormal_asymmetric_DFE_0009987.InferDFE.godambe"
dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_lognormal --dfe-popt "$dfe_popt_file" --cache2d "$cache2d_file" --output "$output_file"


# Lognormal
# Error: None!
# dfe_popt_file="mmd_FRA_mmd_IRA_1d_lognormal_DFE_0009987.InferDFE.bestfits"
# output_file="mmd_FRA_mmd_IRA_1d_lognormal_DFE_0009987.InferDFE.godambe"
# dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d lognormal --dfe-popt "$dfe_popt_file" --cache1d "$cache1d_file" --output "$output_file"


# Symmetric Bivariate Independent Gamma
# Error: None!
# dfe_popt_file="mmd_FRA_mmd_IRA_2d_biv_ind_gamma_symmetric_DFE_0009987.InferDFE.bestfits"
# output_file="mmd_FRA_mmd_IRA_2d_biv_ind_gamma_symmetric_DFE_0009987.InferDFE.godambe"
# dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_ind_gamma --dfe-popt "$dfe_popt_file" --cache2d "$cache2d_file" --output "$output_file"


# Asymmetric Bivariate Independent Gamma
# Error: None!
# dfe_popt_file="mmd_FRA_mmd_IRA_2d_biv_ind_gamma_asymmetric_DFE_0009987.InferDFE.bestfits"
# output_file="mmd_FRA_mmd_IRA_2d_biv_ind_gamma_asymmetric_DFE_0009987.InferDFE.godambe"
# dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_ind_gamma --dfe-popt "$dfe_popt_file" --cache2d "$cache2d_file" --output "$output_file"

# Gamma
# Error: None!
# dfe_popt_file="mmd_FRA_mmd_IRA_2d_gamma_DFE_0009987.InferDFE.bestfits"
# output_file="mmd_FRA_mmd_IRA_2d_gamma_DFE_0009987.InferDFE.godambe"
# dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d gamma --dfe-popt "$dfe_popt_file" --cache1d "$cache1d_file" --output "$output_file"
