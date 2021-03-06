
all:
	#gcc -c test_c_function.c weak_strong_4d_c.c transverse_field_gauss_round.c
	#gcc -c -std=c99 Faddeeva.c 
	#gcc -std=c99 -o test_c_function test_c_function.o weak_strong_4d_c.o Faddeeva.o transverse_field_gauss_round.o -lm
	f2py -m boost_sixtrack -c boost_sixtrack.f
	f2py -m beambeam_force_sixtrack -c beambeam_force_sixtrack.f 
	f2py -m sigmas_sixtrak -c sigmas_sixtrak.f
	f2py -m slicing_sixtrack -c slicing_sixtrack.f
	f2py -m full_interaction_sixtrack -c full_interaction_sixtrack.f 
	gcc -shared -lm -fPIC -Ifrom_sixtracklib -o cFromsixtracklib.so cFromSixtracklib.c
	gcc -shared -lm -fPIC -Icsrc -Ifrom_sixtracklib -o cBB6D.so cBB6D.c
 


