GALACTIC_CENTRE="4.6498509244036468,-0.50628171572274738"
CRAB_PULSAR="1.4596750675891823,0.38422482944113812"
rm -f *.dat *.png
time ../bin/radiometer --injection-ra-dec ${GALACTIC_CENTRE} --injection-ra-dec ${CRAB_PULSAR}
#../bin/plotsky --injection-ra-dec ${GALACTIC_CENTRE} --injection-ra-dec ${CRAB_PULSAR} power_td_*.dat
