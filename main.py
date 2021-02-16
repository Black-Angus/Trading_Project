from SABR_function import *
from input import*

time, tenor = maturite(dfvols)
strike = swaprates
MKT = vol(dfvols)
spot = swaprates

alpha, beta, rho, nu = calibration(starting_guess, f, strike, time, MKT)
matrix_vol = SABR_vol_matrix(alpha, beta, rho, nu, f, strike, time)
prices_call = bachelier_matrix(matrix_vol, spot, strike, f, time, "call")
prices_put = bachelier_matrix(matrix_vol, spot, strike, f, time, "put")
