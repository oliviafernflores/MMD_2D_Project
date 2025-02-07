#!/bin/bash
#SBATCH --job-name=GO_DFE_stats
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=10
#SBATCH --ntasks=50
#SBATCH --time=24:00:00


# Base directory
# base_dir="/groups/rgutenk/oliviafernflores/gene_ontology/GO_DFE_Inference"
base_dir="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference"

# Loop through each directory in the base directory
for dir in "$base_dir"/*; do
    if [ -d "$dir" ]; then  # Check if it's a directory
        # Extract the number from the directory name
        dir_name="${dir##*/}"  # Get just the directory name
        number="${dir_name#*_}"  # Remove everything up to the first underscore
        number="${number%_*}"  # Remove everything from the last underscore onward


        # Define the input and output files based on the directory name
        fs_file="$dir/${number}_dict_Mmd_FRA_Mmd_IRA_SFS_nsyn_unfolded.fs"
        
        # Additional check for directories with a different naming structure
        if [[ "$dir_name" == "GO_with_missing_data" ]]; then
            fs_file="$dir/data_dict_Mmd_FRA_Mmd_IRA_SFS_nsyn_unfolded.fs"
        fi

        cache1d_file="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/FRA_IRA/mmd_FRA_mmd_IRA_1d_cache.bpkl"
        cache2d_file="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/FRA_IRA/mmd_FRA_mmd_IRA_2d_cache.bpkl"
        boots_syn_dir="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_demography/IRA_FRA/bootstraps"
        boots_nsyn_dir="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/GO_bootstraps/GO_{number}_bootstraps"

        # cache1d_file="/groups/rgutenk/oliviafernflores/IRA_FRA/mmd_FRA_mmd_IRA_1d_cache.bpkl"
        # cache2d_file="/groups/rgutenk/oliviafernflores/IRA_FRA/mmd_FRA_mmd_IRA_2d_cache.bpkl"
        # boots_syn_dir="/groups/rgutenk/oliviafernflores/IRA_FRA/bootstraps/"
        # boots_nsyn_dir="/groups/rgutenk/oliviafernflores/gene_ontology/GO_bootstraps/GO_{number}_bootstraps/"

        # Check if the fs file exists before running the commands
        if [ -f "$fs_file" ]; then
            # Run the commands with updated output filenames, saving them in the same directory

            # dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_lognormal --dfe-popt "$dir/mmd_FRA_mmd_IRA_2d_biv_lognormal_DFE_${number}.InferDFE.bestfits" --cache2d "$cache2d_file" --output "$dir/mmd_FRA_mmd_IRA_2d_biv_lognormal_DFE_${number}.params.godambe"

            dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d lognormal --dfe-popt "$dir/mmd_FRA_mmd_IRA_1d_lognormal_DFE_${number}.InferDFE.bestfits" --cache1d "$cache1d_file" --output "$dir/mmd_FRA_mmd_IRA_1d_lognormal_DFE_${number}.params.godambe"

            # dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_lognormal --dfe-popt "$dir/mmd_FRA_mmd_IRA_2d_biv_lognormal_asymmetric_DFE_${number}.InferDFE.bestfits" --cache2d "$cache2d_file" --output "$dir/mmd_FRA_mmd_IRA_2d_biv_lognormal_asymmetric_DFE_${number}.params.godambe"

            dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_ind_gamma --dfe-popt "$dir/mmd_FRA_mmd_IRA_2d_biv_ind_gamma_asymmetric_DFE_${number}.InferDFE.bestfits" --cache2d "$cache2d_file" --output "$dir/mmd_FRA_mmd_IRA_2d_biv_ind_gamma_asymmetric_DFE_${number}.params.godambe"

            dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d biv_ind_gamma --dfe-popt "$dir/mmd_FRA_mmd_IRA_2d_biv_ind_gamma_symmetric_DFE_${number}.InferDFE.bestfits" --cache2d "$cache2d_file" --output "$dir/mmd_FRA_mmd_IRA_2d_biv_ind_gamma_symmetric_DFE_${number}.params.godambe"

            dadi-cli StatDFE --bootstrapping-nonsynonymous-dir "$boots_nsyn_dir" --bootstrapping-synonymous-dir "$boots_syn_dir" --fs "$fs_file" --pdf2d gamma --dfe-popt "$dir/mmd_FRA_mmd_IRA_2d_gamma_DFE_${number}.InferDFE.bestfits" --cache1d "$cache1d_file" --output "$dir/mmd_FRA_mmd_IRA_2d_gamma_DFE_${number}.params.godambe"

        else
            echo "Warning: File $fs_file not found. Skipping directory $dir."
        fi
    fi
done
