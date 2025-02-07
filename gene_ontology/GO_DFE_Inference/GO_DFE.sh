#!/bin/bash
#SBATCH --job-name=GO_DFE_inference
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=24:00:00
#SBATCH --array=1

# Base directory
base_dir="/groups/rgutenk/oliviafernflores/gene_ontology/GO_DFE_Inference"

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

        cache1d_file="/groups/rgutenk/oliviafernflores/IRA_FRA/mmd_FRA_mmd_IRA_1d_cache.bpkl"
        cache2d_file="/groups/rgutenk/oliviafernflores/IRA_FRA/mmd_FRA_mmd_IRA_2d_cache.bpkl"
        demo_popt_file="/groups/rgutenk/oliviafernflores/IRA_FRA/demo_results/IRA_FRA_im_pre_inbreeding_demo_fits_combined.txt"

        # Check if the fs file exists before running the commands
        if [ -f "$fs_file" ]; then
            # Run the commands with updated output filenames, saving them in the same directory
            dadi-cli InferDFE --fs "$fs_file" --cache2d "$cache2d_file" --pdf2d biv_lognormal --p0 1 1 .5 .5 --lbounds -10 0.01 0.001 0 --ubounds 100 100 0.999 0.5 --demo-popt "$demo_popt_file" --ratio 2.4 --output "$dir/mmd_FRA_mmd_IRA_2d_biv_lognormal_DFE_$number" --optimizations 15 --maxeval 400 --check-convergence 10

            dadi-cli InferDFE --fs "$fs_file" --cache1d "$cache1d_file" --pdf2d lognormal --p0 1 1 .5 --lbounds -10 0.01 0 --ubounds 100 100 0.5 --demo-popt "$demo_popt_file" --ratio 2.4 --output "$dir/mmd_FRA_mmd_IRA_1d_lognormal_DFE_$number" --optimizations 15 --maxeval 400 --check-convergence 10

            dadi-cli InferDFE --fs "$fs_file" --cache2d "$cache2d_file" --pdf2d biv_lognormal --p0 1 1 1 1 .5 .5 --lbounds -10 -10 0.01 0.01 0.001 0 --ubounds 100 100 100 100 0.999 0.5 --demo-popt "$demo_popt_file" --ratio 2.4 --output "$dir/mmd_FRA_mmd_IRA_2d_biv_lognormal_asymmetric_DFE_$number" --optimizations 15 --maxeval 400 --check-convergence 10

            dadi-cli InferDFE --fs "$fs_file" --cache2d "$cache2d_file" --pdf2d biv_ind_gamma --p0 1 1 1 1 .5 --lbounds 0.010 0.010 0.01 0.01 0 --ubounds 100 1000 100 1000 0.999 --demo-popt "$demo_popt_file" --ratio 2.4 --output "$dir/mmd_FRA_mmd_IRA_2d_biv_ind_gamma_asymmetric_DFE_$number" --optimizations 20 --maxeval 400 --check-convergence 10

            dadi-cli InferDFE --fs "$fs_file" --cache2d "$cache2d_file" --pdf2d biv_ind_gamma --p0 1 1 .5 --lbounds 0.010 0.010  0 --ubounds 100 1000 0.999 --demo-popt "$demo_popt_file" --ratio 2.4 --output "$dir/mmd_FRA_mmd_IRA_2d_biv_ind_gamma_symmetric_DFE_$number" --optimizations 20 --maxeval 400 --check-convergence 10

            dadi-cli InferDFE --fs "$fs_file" --cache1d "$cache1d_file" --pdf2d gamma --p0 1 1 .5 --lbounds 0.010 0.010  0 --ubounds 100 1000 0.999 --demo-popt "$demo_popt_file" --ratio 2.4 --output "$dir/mmd_FRA_mmd_IRA_2d_gamma_DFE_$number" --optimizations 20 --maxeval 400 --check-convergence 10
        else
            echo "Warning: File $fs_file not found. Skipping directory $dir."
        fi
    fi
done
