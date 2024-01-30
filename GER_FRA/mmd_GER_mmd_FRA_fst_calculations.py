import dadi


data_fs = dadi.Spectrum.from_file('FRA_GER_syn_unfolded.fs')
fst = data_fs.Fst()
print(fst)