#! /usr/bin/python
#
#  Script for performing
#  cylindrical matching between
#  observed and mock clusters.
#

import numpy as np

from functions import errors
from pycymatch_lib import *


##
#  Code Main.
def main():

    ##
    # Read arguments
    opts = pycymatch_opts.get_opts().parse_args()
    opts.input_files[1] = '/Users/mahaixia/VIPERS_MOCK/VIPERS_W4_VIPERSLIKE_MOCKS/mock_W4_001_spec_VAC.dat'
    opts.input_files[0] = '/Users/mahaixia/Documents/GitHub/FOF/Number-Rz/VIPERS_W4.txt_clusters_0.15_0.009_spec.dat'
    opts.obs_cols = 11
    opts.mock_cols = 11

    ##
    # Check input files
    for file in opts.input_files:
        errors.file_name_error(file)

    pycymatch_extra.h_line()

    ##
    # Read mock halo catalogue
    mock = pycymatch_io.read_mock(opts)

    ##
    # Read observed catalogue
    obs = pycymatch_io.read_obs(opts)

    pycymatch_extra.h_line()

    ##
    # Find Matches
    matches = pycymatch_match2.find_matches(mock, obs, 2.0 * opts.delta_z,
                                            opts)
    pycymatch_io.print_matches(mock[matches[2]], obs, matches[0], matches[1],
                               opts)

    ##
    # Define completeness and purity of sample
    c_matrix = pycymatch_cp.get_completeness(mock, matches, opts)
    p_matrix = pycymatch_cp.get_purity(obs, matches, opts)

    ##
    # Define mass-observable matrix for matched objects
    matrix, hm_matrix, ranges = pycymatch_match2.mo_matrix(mock, obs, matches,
                                                           opts)

    pycymatch_extra.h_line()

    ##
    # Make plots
    pycymatch_plot.make_plots(matrix, hm_matrix, ranges, opts)
    pycymatch_plot.plot_complete(c_matrix[0], 'mass', opts)
    pycymatch_plot.plot_complete(c_matrix[1], 'ngal', opts)
    pycymatch_plot.plot_pure(p_matrix, opts)

    ##
    # Save matrix to file
    pycymatch_io.print_matrix(matrix, hm_matrix, ranges, opts)

    pycymatch_extra.h_line()

if __name__ == "__main__":
    main()
