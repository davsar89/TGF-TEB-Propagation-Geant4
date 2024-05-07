

import os
import glob
import numpy as np
# import matplotlib.pyplot as plt
import argparse


class SimData:
	ph = {}
	elec = {}
	posi = {}


def read_Geant4_output(filename, particle_type):

	## Open file
	file = open(filename, 'r')
	file_line_by_line = file.readlines()

	## Collect data (line-by-line)
	data_array = []
	for line in file_line_by_line:
		line_as_array = line.strip('\n').split()
		data_array.append(line_as_array)
	file.close()

	## Transpose data (index by columns)
	data_array = [list(i) for i in zip(*data_array)]

	## Organize data in dictionary
	data_dict = _organize_data(data_array, particle_type)

	return data_dict


def _organize_data(data_array, selection):

	particle_map = {'ph':22, 'e-':11, 'e+':-11}

	seed = np.asarray(data_array[0], dtype=int)            # random number seed of the run
	src_lat = np.asarray(data_array[1], dtype=float)       # deg
	src_lon = np.asarray(data_array[2], dtype=float)       # deg
	src_alt = np.asarray(data_array[3], dtype=float)       # km
	beam_dist = np.asarray(data_array[4], dtype=str)       # 1 = Gaussian | 0 = Uniform
	beam_open = np.asarray(data_array[5], dtype=float)     # deg: Gaussian = sigma | Uniform = half-cone
	beam_tilt = np.asarray(data_array[6], dtype=float)     # deg
	num_ph_sampled = np.asarray(data_array[7], dtype=int)  # number of TGF photons sampled when record is made
	creation_type = np.asarray(data_array[8], dtype=int)   # 1 = primary | >1 = secondary
	particle_type = np.asarray(data_array[9], dtype=int)   # 22 = ph | 11 = e- | -11 = e
	time = np.asarray(data_array[10], dtype=float)         # microsec wrt TGF photon sampling time
	energy = np.asarray(data_array[11], dtype=float)       # keV
	lat = np.asarray(data_array[12], dtype=float)          # deg
	lon = np.asarray(data_array[13], dtype=float)          # deg
	alt = np.asarray(data_array[14], dtype=float)          # km
	dist = np.asarray(data_array[15], dtype=float)         # km: distance wrt TGF source position
	pos_x = np.asarray(data_array[16], dtype=float)        # position x-coord
	pos_y = np.asarray(data_array[17], dtype=float)        # position y-coord
	pos_z = np.asarray(data_array[18], dtype=float)        # position z-coord
	mom_x = np.asarray(data_array[19], dtype=float)        # momentum x-coord
	mom_y = np.asarray(data_array[20], dtype=float)        # momentum y-coord
	mom_z = np.asarray(data_array[21], dtype=float)        # momentum z-coord

	beam_dist[beam_dist=='0'] = 'uniform'
	beam_dist[beam_dist=='1'] = 'gaussian'

	index = np.where(particle_type==particle_map[selection])[0]

	data_dict = {}
	data_dict['seed'] = seed[index]
	data_dict['src_lat'] = src_lat[index]
	data_dict['src_lon'] = src_lon[index]
	data_dict['src_alt'] = src_alt[index]
	data_dict['beam_dist'] = beam_dist[index]
	data_dict['beam_open'] = beam_open[index]
	data_dict['beam_tilt'] = beam_tilt[index]
	data_dict['num_ph_sampled'] = num_ph_sampled[index]
	data_dict['creation_type'] = creation_type[index]
	data_dict['particle_type'] = particle_type[index]
	data_dict['time'] = time[index]
	data_dict['energy'] = energy[index]
	data_dict['lat'] = lat[index]
	data_dict['lon'] = lon[index]
	data_dict['alt'] = alt[index]
	data_dict['dist'] = dist[index]
	data_dict['pos_x'] = pos_x[index]
	data_dict['pos_y'] = pos_y[index]
	data_dict['pos_z'] = pos_z[index]
	data_dict['mom_x'] = mom_x[index]
	data_dict['mom_y'] = mom_y[index]
	data_dict['mom_z'] = mom_z[index]

	return data_dict

def write_file(filename, string):
	file = open(filename, 'a+')
	file.write(string)
	file.close()
	return



def run(args):

	## make output directory (if neccessary)
	os.makedirs(args.outDir, exist_ok=True)

	## input file names
	ph_file = glob.glob(os.path.join(args.inDir, 'detPhotons*.out'))[0]
	lep_file = glob.glob(os.path.join(args.inDir, 'detLeptons*.out'))[0]


	## output file names
	elec_posi_ratio_filename = os.path.join(args.outDir, 'elec_posi_ratio.txt')
	elec_ener_mom_filename = os.path.join(args.outDir, 'electrons_ener_momentums.txt')
	posi_ener_mom_filename = os.path.join(args.outDir, 'positrons_ener_momentums.txt')
	ph_ener_spec_filename = os.path.join(args.outDir, 'photons_spec.txt')
	elec_ener_spec_filename = os.path.join(args.outDir, 'electrons_spec.txt')
	posi_ener_spec_filename = os.path.join(args.outDir, 'positrons_spec.txt')

	spec_plot_filename = os.path.join(args.outDir, 'spectra.png')

	## collect data
	data = SimData()
	data.ph = read_Geant4_output(ph_file, 'ph')
	data.elec = read_Geant4_output(lep_file, 'e+')
	data.posi = read_Geant4_output(lep_file, 'e-')

	## electron-positron ratio ##

	elec_posi_ratio = '{:.3e}'.format(len(data.posi['particle_type'])/len(data.elec['particle_type']))
	write_file(elec_posi_ratio_filename, elec_posi_ratio)

	## energy & momentum ##

	elec_ener_mom = ''
	for i in range(data.elec['energy'].size):
		elec_ener_mom += '{} {} {} {}\n'.format(data.elec['energy'][i], data.elec['mom_x'][i], data.elec['mom_y'][i], data.elec['mom_z'][i])
	write_file(elec_ener_mom_filename, elec_ener_mom)

	posi_ener_mom = ''
	for i in range(data.posi['energy'].size):
		posi_ener_mom += '{} {} {} {}\n'.format(data.posi['energy'][i], data.posi['mom_x'][i], data.posi['mom_y'][i], data.posi['mom_z'][i])
	write_file(posi_ener_mom_filename, posi_ener_mom)

	## spectra ##

	## log-space BGO energy bins
	energy_bins = np.logspace(np.log10(200), np.log10(40_000), num=129)
	## returns index value of bin that that value corresponds to [indexing starts at 1]
	ph_index = np.digitize(data.ph['energy'], energy_bins, right=True)
	elec_index = np.digitize(data.elec['energy'], energy_bins, right=True)
	posi_index = np.digitize(data.posi['energy'], energy_bins, right=True)
	## counts per bin 
	ph_cnts_per_bin = np.zeros(128)
	elec_cnts_per_bin = np.zeros(128)
	posi_cnts_per_bin = np.zeros(128)
	for index in ph_index:
		ph_cnts_per_bin[index-1] += 1
	for index in elec_index:
		elec_cnts_per_bin[index-1] += 1
	for index in posi_index:
		posi_cnts_per_bin[index-1] += 1
	## dnde
	ph_dnde = ph_cnts_per_bin/(energy_bins[1:]-energy_bins[:-1])
	elec_dnde = elec_cnts_per_bin/(energy_bins[1:]-energy_bins[:-1])
	posi_dnde = posi_cnts_per_bin/(energy_bins[1:]-energy_bins[:-1])
	## normalize
	ph_dnde = ph_dnde/np.mean(ph_dnde)
	elec_dnde = elec_dnde/np.mean(elec_dnde)
	posi_dnde = posi_dnde/np.mean(posi_dnde)

	ph_ener_spec = ''
	for i in range(ph_dnde.size):
		ph_ener_spec += '{:.3f} {:.5f}\n'.format(energy_bins[:-1][i], ph_dnde[i])
	write_file(ph_ener_spec_filename, ph_ener_spec)

	elec_ener_spec = ''
	for i in range(elec_dnde.size):
		elec_ener_spec += '{:.3f} {:.5f}\n'.format(energy_bins[:-1][i], elec_dnde[i])
	write_file(elec_ener_spec_filename, elec_ener_spec)

	posi_ener_spec = ''
	for i in range(posi_dnde.size):
		posi_ener_spec += '{:.3f} {:.5f}\n'.format(energy_bins[:-1][i], posi_dnde[i])
	write_file(posi_ener_spec_filename, posi_ener_spec)

	## plot
	# plt.step(energy_bins[:-1], ph_dnde, label='ph')
	# plt.step(energy_bins[:-1], elec_dnde, label='e-')
	# plt.step(energy_bins[:-1], posi_dnde, label='e+')
	# plt.xlabel('Energy (keV)')
	# plt.ylabel('DnDe (normalized)')
	# plt.xscale('log')
	# plt.yscale('log')
	# plt.legend()
	# plt.savefig(spec_plot_filename)





if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('tgf_id', type=str)
	args = parser.parse_args()

	args.inDir = '.'
#	args.inDir = './TGF-TEB-Propagation-Geant4_subfolders/build/output_ascii/541_15_30_Gaussian_0'
#	args.outDir = './TEB_MASS_MODEL_RESPONSE_SL_edit/build/input/'+args.tgf_id+'/'
	args.outDir = '.'
	run(args)










