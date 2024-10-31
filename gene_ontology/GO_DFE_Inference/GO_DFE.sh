#!/usr/bin/env bash

# Base directory where the GO directories are located
base_dir="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/GO_use_for_DFE"

# Check if the base directory exists
if [ ! -d "$base_dir" ]; then
    echo "Base directory not found: $base_dir"
    exit 1
fi

# Loop through each directory matching the pattern GO_*
for dir in "$base_dir"/GO_*; do
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"

        # Loop through each .fs file in the directory
        for fs_file in "$dir"/*_nsyn_unfolded.fs; do
            if [ -f "$fs_file" ]; then
                echo "Found fs file: $fs_file"

                # Extract the population pairs from the file name
                populations=$(basename "$fs_file" | sed -E 's/.*_Mmd_([A-Z]+)_Mmd_([A-Z]+)_SFS.*/\1_\2/')
                
                # Check for both orders
                if [[ "$populations" =~ (FRA_GER|GER_FRA) ]]; then
                    demo_popt_path="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_demography/FRA_GER/demo_results/FRA_GER_im_inbreeding_demo_fits_combined.txt"
                elif [[ "$populations" =~ (IRA_FRA|FRA_IRA) ]]; then
                    demo_popt_path="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_demography/IRA_FRA/demo_results/IRA_FRA_im_pre_inbreeding_demo_fits_combined.txt"
                elif [[ "$populations" =~ (GER_HEL|HEL_GER) ]]; then
                    demo_popt_path="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_demography/GER_HEL/demo_results/GER_HEL_im_pre_inbreeding_demo_fits_combined.txt"
                else
                    echo "No matching demo popt file for populations: $populations"
                    continue
                fi

                # Construct the cache file and output file names
                cache_file="mmd_${populations}_1d_cache.bpkl"

                if [[ "$populations" =~ (FRA_GER|GER_FRA) ]]; then
                    cache_path="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/GER_FRA/mmd_GER_mmd_FRA_1d_cache.bpkl"
                elif [[ "$populations" =~ (IRA_FRA|FRA_IRA) ]]; then
                    cache_path="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/FRA_IRA/mmd_FRA_mmd_IRA_1d_cache.bpkl"
                elif [[ "$populations" =~ (GER_HEL|HEL_GER) ]]; then
                    cache_path="/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/HEL_GER/mmd_HEL_mmd_GER_1d_cache.bpkl"                else
                    echo "No matching cache file for populations: $populations"
                    continue
                fi

                output_file="$dir/${populations}_biv_lognormal_DFE_dadi_cli.txt"

                # Check if cache and demo popt files exist
                if [ ! -f "$cache_path" ]; then
                    echo "Cache file not found: $cache_path"
                    continue
                fi

                if [ ! -f "$demo_popt_path" ]; then
                    echo "Demo popt file not found: $demo_popt_path"
                    continue
                fi

                # Run the command
                dadi-cli InferDFE --fs "$fs_file" \
                    --cache2d "$cache_path" \
                    --pdf2d biv_lognormal \
                    --p0 1 1 0.5 0.05 \
                    --lbounds "-10" 0.01 0.001 0 \
                    --ubounds 10 10 0.999 0.5 \
                    --demo-popt "$demo_popt_path" \
                    --ratio 2.4 \
                    --output "$output_file" \
                    --optimizations 15 \
                    --maxeval 400 \
                    --check-convergence 10

                echo "Processed: $fs_file -> $output_file"
            else
                echo "No .fs file found in: $dir"
            fi
        done
    else
        echo "Not a directory: $dir"
    fi
done
