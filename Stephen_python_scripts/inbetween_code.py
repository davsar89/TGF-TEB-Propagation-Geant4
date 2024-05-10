

import os
import glob
import numpy as np
# import matplotlib.pyplot as plt
import argparse
import subprocess
import copy
import sys

#BASE_OUTPUT_FILE_NAME = 'detLeptons_'
BASE_OUTPUT_FILE_NAME = 'detParticles_'

##
idx_RAND_SEED = 0
idx_SOURCE_ALT = 1
idx_SOURCE_OPENING_ANGLE = 2
idx_TILT_ANGLE = 3
idx_event_nb = 4
idx_ID = 5
idx_PDG_NB = 6
idx_time = 7
idx_energy = 8
idx_alt = 9
idx_lat = 10
idx_lon = 11
idx_dist_rad = 12
idx_ecef_x = 13
idx_ecef_y = 14
idx_ecef_z = 15
idx_mom_x = 16
idx_mom_y = 17
idx_mom_z = 18
idx_number_beaming = 19
idx_SOURCE_LAT = 20
idx_SOURCE_LONG = 21
##

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

def get_separator_lat(lat_array):
	# Sort the data
	xx = copy.deepcopy(lat_array)
	lat_sorted = np.sort(xx)
	# Find gaps in the sorted data
	gaps = np.diff(lat_sorted)
	# Find the index of the largest gap
	index = np.argmax(gaps)
	# Calculate the middle x value between the largest gap
	middle_x = (lat_sorted[index] + lat_sorted[index + 1]) / 2
	return middle_x

def _organize_data(data_array, selection):

	particle_map = {'ph':22, 'e-':11, 'e+':-11}

	seed = np.asarray(data_array[idx_RAND_SEED], dtype=int)            # random number seed of the run
	src_lat = np.asarray(data_array[idx_SOURCE_LAT], dtype=float)       # deg
	src_lon = np.asarray(data_array[idx_SOURCE_LONG], dtype=float)       # deg
	src_alt = np.asarray(data_array[idx_SOURCE_ALT], dtype=float)       # km
	beam_dist = np.asarray(data_array[idx_number_beaming], dtype=str)       # 1 = Gaussian | 0 = Uniform
	beam_open = np.asarray(data_array[idx_SOURCE_OPENING_ANGLE], dtype=float)     # deg: Gaussian = sigma | Uniform = half-cone theta
	beam_tilt = np.asarray(data_array[idx_TILT_ANGLE], dtype=float)     # deg
	num_ph_sampled = np.asarray(data_array[idx_event_nb], dtype=int)  # number of TGF photons sampled when record is made
	creation_type = np.asarray(data_array[idx_ID], dtype=int)   # 1 = primary | >1 = secondary
	particle_type = np.asarray(data_array[idx_PDG_NB], dtype=int)   # 22 = ph | 11 = e- | -11 = e
	time = np.asarray(data_array[idx_time], dtype=float)         # microsec wrt TGF photon sampling time
	energy = np.asarray(data_array[idx_energy], dtype=float)       # keV
	lat = np.asarray(data_array[idx_lat], dtype=float)          # deg
	lon = np.asarray(data_array[idx_lon], dtype=float)          # deg
	alt = np.asarray(data_array[idx_alt], dtype=float)          # km
	dist = np.asarray(data_array[idx_dist_rad], dtype=float)         # km: distance wrt TGF source position
	pos_x = np.asarray(data_array[idx_ecef_x], dtype=float)        # position x-coord
	pos_y = np.asarray(data_array[idx_ecef_y], dtype=float)        # position y-coord
	pos_z = np.asarray(data_array[idx_ecef_z], dtype=float)        # position z-coord
	mom_x = np.asarray(data_array[idx_mom_x], dtype=float)        # momentum x-coord
	mom_y = np.asarray(data_array[idx_mom_y], dtype=float)        # momentum y-coord
	mom_z = np.asarray(data_array[idx_mom_z], dtype=float)        # momentum z-coord

	beam_dist[beam_dist=='0'] = 'uniform'
	beam_dist[beam_dist=='1'] = 'gaussian'

	separator_lat = get_separator_lat(lat)

	kept_indices = np.where((particle_type==particle_map[selection]) & (lat>separator_lat))[0]

	data_dict = {}
	data_dict['seed'] = seed[kept_indices]
	data_dict['src_lat'] = src_lat[kept_indices]
	data_dict['src_lon'] = src_lon[kept_indices]
	data_dict['src_alt'] = src_alt[kept_indices]
	data_dict['beam_dist'] = beam_dist[kept_indices]
	data_dict['beam_open'] = beam_open[kept_indices]
	data_dict['beam_tilt'] = beam_tilt[kept_indices]
	data_dict['num_ph_sampled'] = num_ph_sampled[kept_indices]
	data_dict['creation_type'] = creation_type[kept_indices]
	data_dict['particle_type'] = particle_type[kept_indices]
	data_dict['time'] = time[kept_indices]
	data_dict['energy'] = energy[kept_indices]
	data_dict['lat'] = lat[kept_indices]
	data_dict['lon'] = lon[kept_indices]
	data_dict['alt'] = alt[kept_indices]
	data_dict['dist'] = dist[kept_indices]
	data_dict['pos_x'] = pos_x[kept_indices]
	data_dict['pos_y'] = pos_y[kept_indices]
	data_dict['pos_z'] = pos_z[kept_indices]
	data_dict['mom_x'] = mom_x[kept_indices]
	data_dict['mom_y'] = mom_y[kept_indices]
	data_dict['mom_z'] = mom_z[kept_indices]

	return data_dict

def write_file(filename, string):
	file = open(filename, 'w+')
	file.write(string)
	file.close()
	print(f"Finished writing {filename}.")
	return

def run(args):

	## make output directory (if neccessary)
	os.makedirs(args.outDir, exist_ok=True)

	## input file name
	command0 = f'rm -rf {args.inDir}/fused.out'
	subprocess.run(command0, shell=True, check=True)
	command = f'cat {args.inDir}/*.out > {args.inDir}/fused.out'
	subprocess.run(command, shell=True, check=True)

	concatenated_file = glob.glob(os.path.join(args.inDir, f'fused.out'))[0]

	## output file names
	elec_posi_ratio_filename = os.path.join(args.outDir, 'elec_posi_ratio.txt')
	elec_ener_mom_filename = os.path.join(args.outDir, 'electrons_ener_momentums.txt')
	posi_ener_mom_filename = os.path.join(args.outDir, 'positrons_ener_momentums.txt')
	ph_ener_spec_filename = os.path.join(args.outDir, 'photons_spec.txt')
	elec_ener_spec_filename = os.path.join(args.outDir, 'electrons_spec.txt')
	posi_ener_spec_filename = os.path.join(args.outDir, 'positrons_spec.txt')

	## collect data
	data = SimData()
	if not args.ignore_photons:
		data.ph = read_Geant4_output(concatenated_file, 'ph')
	data.elec = read_Geant4_output(concatenated_file, 'e-')
	data.posi = read_Geant4_output(concatenated_file, 'e+')

	## electron-positron ratio ##

	ep_r = len(data.posi['particle_type']) / (len(data.elec['particle_type']) + len(data.posi['particle_type']))
	elec_posi_ratio = "{:.5e}".format(ep_r)
	write_file(elec_posi_ratio_filename, elec_posi_ratio)

	## energy & momentum ##

	elec_ener_mom_str = ''
	for i in range(data.elec['energy'].size):
		elec_ener_mom_str += '{:.5e} {:.5e} {:.5e} {:.5e}\n'.format(data.elec['energy'][i], data.elec['mom_x'][i], data.elec['mom_y'][i], data.elec['mom_z'][i])
	write_file(elec_ener_mom_filename, elec_ener_mom_str)

	posi_ener_mom_str = ''
	for i in range(data.posi['energy'].size):
		posi_ener_mom_str += '{:.5e} {:.5e} {:.5e} {:.5e}\n'.format(data.posi['energy'][i], data.posi['mom_x'][i], data.posi['mom_y'][i], data.posi['mom_z'][i])
	write_file(posi_ener_mom_filename, posi_ener_mom_str)

	## spectra ##

	## log-space BGO energy bins
	energy_bins = np.logspace(np.log10(200), np.log10(40_000), num=129)
	## returns index value of bin that that value corresponds to [indexing starts at 1]
	
	# Calculate histograms
	if not args.ignore_photons:
		ph_hist, bin_edges = np.histogram(data.ph['energy'], bins=energy_bins)
	elec_hist, bin_edges = np.histogram(data.elec['energy'], bins=energy_bins)
	posi_hist, bin_edges = np.histogram(data.posi['energy'], bins=energy_bins)

	# Calculate bin widths
	bin_widths = np.diff(bin_edges)

	# Normalize by the bin size to get a probability density function
	if not args.ignore_photons:
		ph_pdf = ph_hist / (sum(ph_hist) * bin_widths)
	elec_pdf = elec_hist / (sum(elec_hist) * bin_widths)
	posi_pdf = posi_hist / (sum(posi_hist) * bin_widths)

	if not args.ignore_photons:
		ph_ener_spec_str = '{:.5e} {:.5e}\n'.format(energy_bins[0], 0)
		for i in range(ph_pdf.size):
			ph_ener_spec_str += '{:.5e} {:.5e}\n'.format(energy_bins[i+1], ph_pdf[i])
		write_file(ph_ener_spec_filename, ph_ener_spec_str)

	elec_ener_spec_str = '{:.5e} {:.5e}\n'.format(energy_bins[0], 0)
	for i in range(elec_pdf.size):
		elec_ener_spec_str += '{:.5e} {:.5e}\n'.format(energy_bins[i+1], elec_pdf[i])
	write_file(elec_ener_spec_filename, elec_ener_spec_str)

	posi_ener_spec_str = '{:.5e} {:.5e}\n'.format(energy_bins[0], 0)
	for i in range(posi_pdf.size):
		posi_ener_spec_str += '{:.5e} {:.5e}\n'.format(energy_bins[i+1], posi_pdf[i])
	write_file(posi_ener_spec_filename, posi_ener_spec_str)

	## plot spectra
	# spec_plot_filename = os.path.join(args.outDir, 'spectra.png')
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
	parser.add_argument('--tgf_id', type=int, default=99999, help='Optional argument with a default value')
	args = parser.parse_args()

	if args.tgf_id == 99999:
		print(f'Warning: Using the default value for --tgf_id = {args.tgf_id}', file=sys.stderr)
	else:
		print(f'Using --tgf_id = {args.tgf_id}')

	args.inDir = './input/'
#	args.inDir = './TGF-TEB-Propagation-Geant4_subfolders/build/output_ascii/541_15_30_Gaussian_0'
#	args.outDir = './TEB_MASS_MODEL_RESPONSE_SL_edit/build/input/'+args.tgf_id+'/'
	args.outDir = './output/'
	args.ignore_photons = True
	run(args)

	print("Done processing.")










