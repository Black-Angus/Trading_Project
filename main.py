from SABR_function import *
time = maturite(dfvols)

calibration(starting_guess, F, K, expiries, MKT)

SABR_vol_matrix(alpha, beta, rho, nu, F, K, expiries, MKT)

outvol.close()
vol_diff.close()
parameters.close()
