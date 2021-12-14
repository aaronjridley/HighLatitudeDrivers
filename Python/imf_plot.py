#!/usr/bin/env python

import sys
import datetime as dt
import numpy as np
import argparse
from swmf_imf import *

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Plot SWMF IMF files')
    parser.add_argument('files', metavar = 'file', nargs = '+', \
                        help = 'Files to process')
    parser.add_argument('-mach', \
                        help='Plot Alfven Mach number', \
                        action="store_true")

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

args = parse_args()

file = args.files[0]
data = read_swmf_file(file)

plot_imf(data)

