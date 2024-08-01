#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="mmd_flip_calls_dd"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=72:00:00

import pickle, dadi, pandas as pd, numpy as np

dd_flip_syn = pickle.load(open('mus_musculus_domesticus.dd_syn_filtered.bpkl', 'rb'))
dd_flip_nsyn = pickle.load(open('mus_musculus_domesticus.dd_nsyn_filtered.bpkl', 'rb'))


for key, value in dd_flip_syn.items():
    calls_dict = value['calls']
    for call_key, call_value in calls_dict.items():
        calls_dict[call_key] = call_value[1], call_value[0]

for key, value in dd_flip_nsyn.items():
    calls_dict = value['calls']
    for call_key, call_value in calls_dict.items():
        calls_dict[call_key] = call_value[1], call_value[0]

with open('mus_musculus_domesticus.dd_syn_correct.bpkl', 'wb') as pick:
    pickle.dump(dd_flip_syn, pick, 2)

with open('mus_musculus_domesticus.dd_nsyn_correct.bpkl', 'wb') as pick:
    pickle.dump(dd_flip_nsyn, pick, 2)