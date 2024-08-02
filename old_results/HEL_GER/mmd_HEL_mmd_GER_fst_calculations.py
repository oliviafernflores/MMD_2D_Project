import dadi


data_fs = dadi.Spectrum.from_file('GER_HEL_syn_unfolded.fs')
fst = data_fs.Fst()
print(fst)