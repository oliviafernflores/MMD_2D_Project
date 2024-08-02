import dadi


data_fs = dadi.Spectrum.from_file('IRA_FRA_syn_unfolded.fs')
fst = data_fs.Fst()
print(fst)