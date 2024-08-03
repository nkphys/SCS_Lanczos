# SCS_Lanczos*
1) This code is written by Nitin Kaushal(kaushalnitin002@gmail.com)
This solves Strongly correlated systems
using the Lanczos Algorithm.
This code also does ED(without Lanczos) for "small enough" systems.

*SCS  stands for Strongly correlated systems.

NOTES for input script:

1) if "Get_overlap_with_basis = true" is given in input, it makes
"Get_Full_Spectrum = true" and "need_few_eig_vecs = true" as well.
2) if "need_few_eig_vecs = true" is given in input, in makes ""Get_Full_Spectrum = true"  as well.
3) "States_to_look = 10 0 y r d r.....", it makes "few_ = 10" as well.
4) remember in "States_to_look = 10 0 y r d r.....", the first state should be always zero.
5) DUSE_COMPLEX flag should be used to do calculations in complex space.
