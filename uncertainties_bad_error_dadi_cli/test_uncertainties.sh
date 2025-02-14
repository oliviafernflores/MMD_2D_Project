boots_nsyn_dir="/Users/olivia/Documents/GO_0009987_nsyn_bootstraps"
boots_syn_dir="/Users/olivia/Documents/IRA_FRA_syn_bootstraps"
fs_file="0009987_Mmd_FRA_Mmd_IRA_SFS_nsyn_unfolded.fs"
cache2d_file="mmd_FRA_mmd_IRA_2d_cache.bpkl"
cache1d_file="mmd_FRA_mmd_IRA_1d_cache.bpkl"


# Lognormal
# Error: None, should be because paths to bootstraps do not exist. Incorrect path for fs or cache(s) does raise error.
dfe_popt_file="mmd_FRA_mmd_IRA_1d_lognormal_DFE_0009987.InferDFE.bestfits"
output_file="mmd_FRA_mmd_IRA_1d_lognormal_DFE_0009987.InferDFE.godambe"
dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d lognormal --dfe-popt "$dfe_popt_file" --cache1d "$cache1d_file" --output "$output_file"

