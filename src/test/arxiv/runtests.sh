# how to generate variants
#
# At 512 Hz:
#
# src/lib/projection.c:
#	projection_matrix_n_elements():
#		N_T = 7 --> buffer = 0
#		N_T = 127 --> buffer = 60
#	projection_matrix_l_max():
#		l_T = 9 --> return l_max
#		l_T = 13 --> return l_max + 4
# src/lib/correlator.c:
#	correlator_power_l_max():
#		l_xi = 17 --> return l_max
#		l_xi = 19 --> return l_max + 2
#
# Similarly for 4096 Hz
#

for d in "127_13_19" "163_69_131" "7_9_17" ; do
	pushd $d

	#rm -f tests_*.dat tests_*.{png,eps} tests.log
	#../../../bin/tests --trials 1 --sample-frequency 512 >tests.log 2>&1
	#../../../bin/tests --trials 1 --sample-frequency 4096 >>tests.log 2>&1
	../../../bin/plottests tests_exact_*.dat

	popd
done
