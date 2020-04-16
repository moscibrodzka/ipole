#!/bin/sh
#
#usage: ./ipole theta freq filename Munit trat_j trat_d sigma_cut
#
# for model_harm3d.c model_harm2d.c model_iharm3d.c etc
./ipole 20. 230e9 HARM3D.001500.h5 1e19 3. 3. 2.

# iharm3d.c
#./ipole 17. 230e9 ../polarized_code_comparizon/mhdmodels/dump_00001096_MADa-0.5.h5 2.87962e+25 1. 160. 1.

#
# Novikov-Thorne disk test freq=1keV
# usage: ./ipole theta freq BHspin
#./ipole 75. 2.41e17 0.99

