#!/usr/bin/env python

import platform
from subprocess import call
from mpi4py import MPI
from functools import partial
from multiprocessing.dummy import Pool
import random
import datetime
import time

# Initialization of MPI routines
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
##

if rank == 0:
    print('Number of threads ' + str(nprocs))

# number of CPU to use (for local run only, not cluster)
NB_CPU_LOCAL = 3

# number of times the list of command is copied (usually usefull for large number of CPU, like for clusters)
NB_super_copy = 2000

nr_to_shoot_per_run = 100000000

# lists of parameters to test
#initial_alt_list = [18.,16.,14.,12.,10.]
#weights_alt = [1.,2.,3.,7.,10.]
initial_alt_list = [15.]
weights_alt = [1.]
#initial_alt_list = [18., 16.]
#beaming_angle_list = [ 30., 40.,50.,60.]
beaming_angle_list = [40.]
weights_ang = [1.]
#beaming_angle_list = [5., 10.,20.,30.,40.,50.,60.]
#weights_ang = [1.,2.,3.,4.,6.,8.,10.]
# tilt_angle_list = [5.,10.,20.,30.,40.,50.,60.]
tilt_angle_list = [0.]
# beaming_type_list=["Gaussian","Uniform"]
beaming_type_list = ["Gaussian"]
source_sigma_time_list = [0.0]

# defining a list of commands to be run in parallel

commands = []

excecutable = './TGF_Propa'

total_weight = 0.0

for i_alt, _ in enumerate(initial_alt_list):
    for i_ang, _ in enumerate(beaming_angle_list):
        weight = float(weights_alt[i_alt]*weights_ang[i_ang])
        total_weight = total_weight + weight

# loops over required initial parameters list

for __ in range(NB_super_copy):
    for i_alt, ini_alt in enumerate(initial_alt_list):
        for i_ang, beaming_ang in enumerate(beaming_angle_list):
            for source_sigmat in source_sigma_time_list:
                for beaming_type in beaming_type_list:
                    for tilt_ang in tilt_angle_list:
                        factor = float(
                            weights_alt[i_alt]*weights_ang[i_ang])/total_weight
                        nb_copies_command = int(round(float(nprocs)*factor))

                        if (nb_copies_command < 1):
                            nb_copies_command = 1

                        if (rank == 0):
                            print(nb_copies_command)

                        for _ in range(nb_copies_command):
                            commands.append(excecutable + ' '+str(nr_to_shoot_per_run)+' '+str(ini_alt)
                                            + ' '+str(beaming_ang)+' ' + str(tilt_ang) + ' '+str(beaming_type)+' '+str(source_sigmat))


# if number of commands is less than NB_CPU_LOCAL, fill the list with extra commands that are duplicate of previous ones
jj = 0
while len(commands) < nprocs:
    commands.append(commands[jj])
    jj += 1
    if jj == len(commands):
        jj = 0

# print(commands)
command_number = len(commands)
if rank == 0:
    print('Number of commands required ' + str(command_number))

time.sleep(6)

####################################
computer_name = platform.node()

if "iftrom" in computer_name:  # run on local (personal) computer

    nb_thread = NB_CPU_LOCAL  # number of threads (cpu) to run

    # Making an array where each element is the list of command for a given thread

    pool = Pool(nb_thread)  # to be always set to 1 for this MPI case
    for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
        if returncode != 0:
            print("%d command failed: %d" % (i, returncode))

else:  # run on computer cluster using MPI

    listy = [[] for i in range(0, nprocs, 1)]

    i = 0
    for j in range(0, command_number, 1):
        listy[i].append(commands[j])
        i = i+1
        if i == nprocs:
            i = 0

    pool = Pool(1)  # to be always set to 1 for this MPI case
    for i, returncode in enumerate(pool.imap(partial(call, shell=True), listy[rank])):
        if returncode != 0:
            print("%d command failed: %d" % (i, returncode))
