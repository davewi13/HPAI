# HPAI
This repository contains the code required to perform statistical inference for epidemiological parameters of HPAI in seabird colonies. The code is directly adapted from that previously applied to African swine fever virus on pig farms. As such there are a few holdovers from that code (e.g. the code references herds rather than colonies as would be more appropriate for seabirds) that could be tidied up if desired but that will not impact the running or usability of the code. The code is set up so that chains can be run in parallel on a computing cluster via slurm. The code was compiled using the GNU C++ compiler as follows

g++ HPAI_inference.cc -O3 -std=c++11 -fmax-errors=3 -o run_model.out
sbatch run_model.sh

An example shell script is included in this repository.
