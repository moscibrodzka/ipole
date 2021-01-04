#!/bin/sh
#
# to run different models makefile needs to substitute model_harm3d.c to your desired model
# for model_harm3d.c model_analytic.c model_iharm3d.c etc
#usage: ./ipole theta freq filename Munit trat_j trat_d sigma_cut
# example
./ipole 20. 230e9 HARM3D.001500.h5 1e19 3. 3. 2.

# when code configured to run the Novikov-Thorne disk test 
# usage: ./ipole theta freq BHspin
# example
#./ipole 75. 2.41e17 0.99
