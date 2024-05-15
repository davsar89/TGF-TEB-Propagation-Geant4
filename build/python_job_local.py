import platform
import subprocess as sp
from subprocess import call
import time
from os import system
from os import walk
from os import listdir
from functools import partial
from multiprocessing.dummy import Pool
import sys

computer_name = platform.node()

NB_CPU = 4
nb_run = 10

################################################################################
# Defining commands to run

# defining the commands to be run in parallel
commands = []
excecutable = './TGF_Propa'

#input parameters
NB = 10000000

SOURCE_ALT = 15
SOURCE_LAT = -7
SOURCE_LONG = 55
SOURCE_SIGMA_TIME = 0
SOURCE_OPENING_ANGLE = 30
TILT_ANGLE = 0
BEAMING_TYPE = 'Gaussian'
RECORD_ALTITUDE = 408

SPECTRUM_MODEL = 0
# 0 for classical RREA 1/Exp(-E/7300)
# 1 for Bowers 2018 reverse positron beam TGF
# 2 for leader Celestin 2015, 60 MV
# 3 for leader Celestin 2015, 160 MV

# DEFINE LIST OF COMMANDS TO RUN ON NB_CPU THREADS
for _ in range(nb_run):
    commands.append(f"{excecutable} {NB} {SOURCE_ALT} {SOURCE_LAT} {SOURCE_LONG} {SOURCE_SIGMA_TIME} {SOURCE_OPENING_ANGLE} {TILT_ANGLE} {BEAMING_TYPE} {RECORD_ALTITUDE} {SPECTRUM_MODEL}")

################################################################################
# LOCAL RUN (uses python multiprocessing library)

nb_thread = NB_CPU  # number of threads (cpu) to use

# Making an array where each element is the list of command for a given thread

command_number = len(commands)

print('Number of commands required ' + str(command_number))

pool = Pool(nb_thread)
for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
    if returncode != 0:
        print("%d command failed: %d" % (i, returncode))
            
