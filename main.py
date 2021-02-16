from SABR_function import *
from input_output import*

time, tenor = maturite(dfvols)
strike = (-200,-150,-100,-50,0,50,100,150,200)
MKT = vol(dfvols)
spot = swaprates
option=[]


alpha, beta, rho, nu = calibration(starting_guess, f, strike, time, MKT)
matrix_vol = SABR_vol_matrix(alpha, beta, rho, nu, f, strike, time)
prices = bachelier_matrix(matrix_vol, spot, strike, f, time, option)