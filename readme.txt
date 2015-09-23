This is an implementation of the Gregory-Loredo algorithm to determine the presence of a periodic signal from a list of arrival times

The implementation consists of 3 functions and 2 scripts to test the algorithm

GL_algorithm.py

This computes the odds-ratio for the presence of a periodic signal, the probability, the most likely spectrum and bin number for details see :

Gregory, P. C. and Thomas. J. Loredo, 1992, "A New Method For The Detection Of A Periodic Signal Of Unknown Shape And Period" in The Astrophysical Journal, Astrophysical J., 398, p.146

compute_bin.py

this computes the bin historgram for a given list of arrival times, a frequency, a phase and a number of bins

simulate_arrival_times.py

simulates periodic rate or constant rate arrival times

test_GL_algorithm.py

This tests the function GL_algorithm with simulated data

test_GL_real_data.py

this uses 2 real data sets (source:http://astrostatistics.psu.edu/datasets/Chandra_flares.html) and does a benchmark on parallel vs. serial execution of the GL_algorithm function

The script looks for files in a folder data in the same directory as the script