import xlrd
import math
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize


def formula(alpha, beta, nu, rho, v, time, logFK):
    A = 1 + (((1 - beta) ** 2 * alpha ** 2) / (24. * (v ** 2)) +
             (alpha * beta * nu * rho) / (4. * v) + ((nu ** 2) * (2 - 3 * (rho ** 2)) / 24.)) * time
    B = 1 + (1 / 24.) * (((1 - beta) * logFK) ** 2) + (1 / 1920.) * (((1 - beta) * logFK) ** 4)
    return A, B


def formula_vect(par, v, time, logFK):
    A = 1 + (((1 - par[1]) ** 2 * par[0] ** 2) / (24. * (v ** 2)) + (par[0] * par[1] * par[3] * par[2]) / (
            4. * v) + ((par[3] ** 2) * (2 - 3 * (par[2] ** 2)) / 24.)) * time
    B = 1 + (1 / 24.) * (((1 - par[1]) * logFK) ** 2) + (1 / 1920.) * (((1 - par[1]) * logFK) ** 4)
    return A, B


def SABR(alpha, beta, rho, nu, forward, strike, time, market_vol):
    logFK = math.log(forward / strike)
    v = (forward * strike) ** ((1 - beta) / 2.)
    a, b = formula(alpha, beta, nu, rho, v, time, logFK)
    if strike <= 0:  # the strike become negative so the smile need to be shifted.
        vol = 0
        diff = 0
    elif forward == strike:  # ATM formula
        vol = (alpha / v) * a
        diff = vol - market_vol
    elif forward != strike:  # not-ATM formula
        z = (nu / alpha) * v * logFK
        x = math.log((math.sqrt(1 - 2 * rho * z + z ** 2) + z - rho) / (1 - rho))
        vol = (nu * logFK * a) / (x * b)
        diff = vol - market_vol

    print(round(vol, 4), '\t', )
    outvol.write('%r;' % round(vol, 4))
    if market_vol == 0:
        diff = 0
        vol_diff.write('%s;' % 'No market data')
    else:
        vol_diff.write('%r;' % round(diff, 4))


def smile(alpha, beta, rho, nu, forward, strike, time, market_vol, i):
    # forward, time and the parameters are scalars,
    # strike and market_vol are vectors,
    # i is the index for tenor/expiry label

    print(label_ten[i], '\t', label_exp[i], '\t')
    outvol.write('%s;%s;' % (label_ten[i], label_exp[i]))
    vol_diff.write('%s;%s;' % (label_ten[i], label_exp[i]))
    parameters.write('%s;%s;' % (label_ten[i], label_exp[i]))

    for j in range(len(strike)):
        if strike[0] <= 0:
            shift(forward, strike)
        SABR(alpha, beta, rho, nu, forward, strike[j], time, market_vol[j])

    print(' ')
    outvol.write('\n')
    vol_diff.write('\n')
    parameters.write('%f;%f;%f;%f;' % (alpha, beta, rho, nu))
    parameters.write('\n')


def SABR_vol_matrix(alpha, beta, rho, nu, forward, strike, time, market_vol):
    # forward, time and the parameters are vectors,
    # strike and market_vol are matrices

    print(' ')
    print((2 + ((num_strikes - 1) / 2)), '       ', 'SABR VOLATILITIES')
    print('  ', '\t', 'strikes:')
    for i in range(num_strikes):
        print(label_strikes[i], '\t')
    print(' ')
    outvol.write('%s;' % 'SABR VOLATILITIES')
    outvol.write('\n')
    vol_diff.write('%s;' % 'VOLATILITY DIFFERENCES')
    vol_diff.write('\n')
    parameters.write('%s;' % 'PARAMETERS')
    parameters.write('\n')
    outvol.write('%s;%s;' % (' ', 'strikes:'))
    vol_diff.write('%s;%s;' % (' ', 'strikes:'))
    for j in range(len(strike_spreads)):
        outvol.write('%s;' % label_strikes[j])
        vol_diff.write('%s;' % label_strikes[j])
    outvol.write('\n')
    vol_diff.write('\n')
    print('tenor', '\t', 'expiry')
    parameters.write('%s;%s;%s;%s;%s;%s' % ('tenor', 'expiry', 'alpha', 'beta', 'rho', 'nu'))
    parameters.write('\n')

    for i in range(len(forward)):
        smile(alpha[i], beta[i], rho[i], nu[i], forward[i], strike[i], time[i], market_vol[i], i)


def shift(forward, strike):
    shift = 0.001 - strike[0]
    for j in range(len(strike)):
        strike[j] = strike[j] + shift
        forward = forward + shift


def objective(par, forward, strike, time, market_vol):
    sum_sq_diff = 0
    if strike[0] <= 0:
        shift(forward, strike)
    for j in range(len(strike)):
        v = (forward * strike[j]) ** ((1 - par[1]) / 2.)
        logFK = math.log(forward / strike[j])
        a, b = formula_vect(par, v, time, logFK)
        if market_vol[j] == 0:
            diff = 0
        elif forward == strike[j]:
            vol = (par[0] / v) * a
            diff = vol - market_vol[j]
        elif forward != strike[j]:
            z = (par[3] / par[0]) * v * logFK
            x = math.log((math.sqrt(1 - 2 * par[2] * z + z ** 2) + z - par[2]) / (1 - par[2]))
            vol = (par[3] * logFK * a) / (x * b)
            diff = vol - market_vol[j]
        sum_sq_diff = sum_sq_diff + diff ** 2
        obj = math.sqrt(sum_sq_diff)
    return obj


def calibration(starting_par, forward, strike, time, market_vol):
    for i in range(len(forward)):
        x0 = starting_par
        bnds = ((0.001, None), (0, 1), (-0.999, 0.999), (0.001, None))
        res = minimize(objective, x0, (forward[i], strike[i], time[i], market_vol[i]), bounds=bnds, method='SLSQP')
        # for a constrained minimization of multivariate scalar functions
        alpha[i] = res.x[0]
        beta[i] = res.x[1]
        rho[i] = res.x[2]
        nu[i] = res.x[3]


def black_scholes(forward, strike, spot, time, market_vol, callorput):
    d1 = (np.log(spot / K) + (forward + 0.5 * market_vol ** 2) * time) / (market_vol * np.sqrt(time))
    d2 = d1() - market_vol * np.sqrt(time)

    if callorput == "call":
        return spot * norm.cdf(d1()) - K * np.exp(-forward * time) * \
               norm.cdf(d2())

    elif calorput == "put":
        return -spot * norm.cdf(-d1()) + strike * np.exp(-forward * time) * \
               norm.cdf(-d2())


def black_scholes_matrix(forward, strike, spot, time, market_vol, callorput):
    print(' ')
    print((2 + ((num_strikes - 1) / 2)), '       ', 'SABR Price')
    print('  ', '\t', 'Price:')
    for i in range(num_strikes):
        print(label_strikes[i], '\t')
    print(' ')
    price.write('%s;' % 'Price')
    price.write('\n')
    price.write('%s;%s;' % (' ', 'strikes:'))
    for j in range(len(strike_spreads)):
        price.write('%s;' % label_strikes[j])
    price.write('\n')
    print('tenor', '\t', 'expiry')
    for i in range(len(forward)):
        black_scholes(forward[i], strike[i],spot[i], time[i], market_vol[i], callorPut[i])


################################################################################################################################################


######## inputs and outputs #########################################

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

######## set swaptions characteristics ###############################

strike_spreads = []
j = 0
while True:
    try:
        strike_spreads.append(int(Market_data.cell(1, 3 + j).value))
        j = j + 1
    except:
        break
num_strikes = len(strike_spreads)

expiries = []
i = 0
while True:
    try:
        expiries.append(Market_data.cell(2 + i, 1).value)
        i = i + 1
    except:
        break

tenors = []
i = 0
while True:
    try:
        tenors.append(Market_data.cell(2 + i, 0).value)
        i = i + 1
    except:
        break

# to create the ATM forward rates
F = []
i = 0
while True:
    try:
        F.append(Market_data.cell(2 + i, 2).value)
        i = i + 1
    except:
        break

# to create the strike grid
K = np.zeros((len(F), num_strikes))
for i in range(len(F)):
    for j in range(num_strikes):
        K[i][j] = F[i] + 0.0001 * (strike_spreads[j])

    # to create market volatilities
MKT = np.zeros((len(F), num_strikes))
for i in range(len(F)):
    for j in range(num_strikes):
        MKT[i][j] = Market_data.cell(2 + i, 3 + j).value

# set starting parameters
starting_guess = np.array([0.001, 0.5, 0, 0.001])
alpha = len(F) * [starting_guess[0]]
beta = len(F) * [starting_guess[1]]
rho = len(F) * [starting_guess[2]]
nu = len(F) * [starting_guess[3]]

######## set labels ###################################################

exp_dates = len(expiries) * [0]
for i in range(len(expiries)):
    if expiries[i] < 1:
        exp_dates[i] = str(int(round(12 * expiries[i]))) + 'm'
    else:
        exp_dates[i] = str(int(round(expiries[i]))) + 'y'
        if expiries[i] - round(expiries[i]) > 0:
            exp_dates[i] = exp_dates[i] + str(int(round((12 * (round(expiries[i], 2) - int(expiries[i])))))) + 'm'
        elif expiries[i] - round(expiries[i]) < 0:
            exp_dates[i] = str(int(round(tenors[i])) - 1) + 'y'
            exp_dates[i] = exp_dates[i] + str(int(round((12 * (round(expiries[i], 2) - int(expiries[i])))))) + 'm'

ten_dates = len(tenors) * [0]
for i in range(len(tenors)):
    if tenors[i] < 1:
        ten_dates[i] = str(int(round(12 * tenors[i]))) + 'm'
    else:
        ten_dates[i] = str(int(round(tenors[i]))) + 'y'
        if tenors[i] - round(tenors[i]) > 0:
            ten_dates[i] = ten_dates[i] + str(int(round((12 * (round(tenors[i], 2) - int(tenors[i])))))) + 'm'
        elif tenors[i] - round(tenors[i]) < 0:
            ten_dates[i] = str(int(round(tenors[i])) - 1) + 'y'
            ten_dates[i] = ten_dates[i] + str(int(round((12 * (round(tenors[i], 2) - int(tenors[i])))))) + 'm'

label_exp = exp_dates
label_ten = ten_dates
label_strikes = num_strikes * [0]
for i in range(num_strikes):
    if strike_spreads[i] == 0:
        label_strikes[i] = 'ATM'
    else:
        label_strikes[i] = str(strike_spreads[i])

######## Call the functions #################################

calibration(starting_guess, F, K, expiries, MKT)

SABR_vol_matrix(alpha, beta, rho, nu, F, K, expiries, MKT)

######## Close output files #################################

outvol.close()
vol_diff.close()
parameters.close()
