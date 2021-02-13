import xlrd
outvol = open('outvol.csv', 'w')  # file output of volatilities
vol_diff = open('vol differences.csv', 'w')  # file output differences between SABR and Market volatilities
parameters = open('parameters.csv', 'w')  # file output parameters
price = open('price.csv', 'w')  # file output price

while True:
    try:
        file_input = xlrd.open_workbook('market_data.xlsx')  # load market data
    except:
        print('Input file is not in the directory!')
    break
Market_data = file_input.sheet_by_name('Swaptions data')  # file input forward rates
