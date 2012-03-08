Limits for RS gravitons

main engine is the roostats_cl95s.C

interface is through simple c++ Limits_tables_Sept15th_CLs.C

cl95 calculator takes the following input:

 Double_t             limit = roostats_cl95(ilum, slum, eff, seff, bck, sbck, n, gauss = false, nuisanceModel, method, 
plotFileName, seed); 
 LimitResult expected_limit = roostats_clm(ilum, slum, eff, seff, bck, sbck, ntoys, nuisanceModel, method, seed); 

i.e, the lumi, signal, bkgd, and nobs, with uncertainties.


for each M1,k/Mpl point, these values need to be calculated, and input to the arrays. 

Limits_tables_Sept15th_CLs.C simply takes these arrays and passes to roostats_cl95s.C


