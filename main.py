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

wb = openpyxl.load_workbook('queue.xls')
sheet = wb.get_sheet_by_name('Sheet1')
sheet.cell(row=(4),column=3).value = alpha
sheet.cell(row=(5),column=3).value = beta
sheet.cell(row=(6),column=3).value = rho
sheet.cell(row=(7),column=3).value = nu
sheet.cell(row=(8),column=3).value = prices_call
sheet.cell(row=(9),column=3).value = prices_put
sheet.cell(row=(10),column=3).value = time
sheet.cell(row=(11),column=3).value = tenor
sheet.cell(row=(12),column=3).value = strike
sheet.cell(row=(13),column=3).value = MKT
sheet.cell(row=(14),column=3).value = spot
wb.save('queue.xls')
