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

wb = openpyxl.load_workbook('Option_3_DataSet.xls')
sheet = wb.get_sheet_by_name('Option_3_DataSet')
for data in enumerate(array):
    worksheet.write_column(2, 9, prices_call)
for data in enumerate(array):
    worksheet.write_column(2, 10, prices_put)
for data in enumerate(array):
    worksheet.write_column(2, 11, time)
for data in enumerate(array):
    worksheet.write_column(2, 12, tenor)
for data in enumerate(array):
    worksheet.write_column(2, 13, strike)
for data in enumerate(array):
    worksheet.write_column(2, 14, MKT)
for data in enumerate(array):
    worksheet.write_column(2, 15, spot)
    
sheet.cell(row=(2),column=5).value = alpha
sheet.cell(row=(2),column=6).value = beta
sheet.cell(row=(2),column=7).value = rho
sheet.cell(row=(2),column=8).value = nu

wb.save('Option_3_DataSet.xls')
