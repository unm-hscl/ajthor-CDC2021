#!/bin/bash

# This script is the entrypoint for the CodeOcean capsule. It calls the
# 'run_all.m' Matlab script. Run this script to generate all figures and tables.
# Results are stored in the 'results' folder.

# Run the MATLAB script
matlab -nodisplay -nosoftwareopengl -r "run_all;"
