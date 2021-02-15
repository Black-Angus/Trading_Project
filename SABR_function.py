from set_label import *
import math
from scipy.stats import norm
from scipy.optimize import minimize

def formula(alpha, beta, nu, rho, v, time, logFK):
    A = 1 + (((1 - beta) ** 2 * alpha ** 2) / (24. * (v ** 2)) +
             (alpha * beta * nu * rho) / (4. * v) + ((nu ** 2) * (2 - 3 * (rho ** 2)) / 24.)) * time
    B = 1 + (1 / 24.) * (((1 - beta) * logFK) ** 2) + (1 / 1920.) * (((1 - beta) * logFK) ** 4)
    return A, B


def SABR(alpha: float, beta: float, rho: float, nu: float, forward: float , strike: float, time: float, market_vol: int):

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


def smile(alpha: float, beta: float, rho: float, nu: float, forward: float , strike: list, time: float, market_vol: list, i:int):

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


def SABR_vol_matrix(alpha: list, beta: list, rho: list, nu: list, forward: list, strike: bytearray, time: list, market_vol: bytearray):

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
        logFK = math.log(forward / strike[j])
        v = (forward * strike[j]) ** ((1 - par[1]) / 2.)
        a, b = formula(par[0], par[1], par[3], par[2], v, time, logFK)
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


def bachelier(market_vol, spot, strike, forward, time, option):

    if option == 'call':
        d = (spot * np.exp(forward * time) - strike) / np.sqrt(market_vol ** 2 / (2 * forward) * (np.exp(2 * forward * time) - 1))

        p = np.exp(-forward * time) * (spot * np.exp(forward * time) - strike) * si.norm.cdf(d) + \
            np.exp(-forward * time) * np.sqrt(market_vol ** 2 / (2 * forward) * (np.exp(2 * forward * time) - 1)) * si.norm.pdf(d)
    if option == 'put':
        p = bachelier(market_vol, spot, strike, forward, time, 'call') - spot + np.exp(-forward * time) * strike

    return p


def bachelier_matrix(forward, strike, spot, time, market_vol, option):
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
        bachelier(forward[i], strike[i], spot[i], time[i], market_vol[i], option[i])


def d1(S, K, T, r, sigma):
    return (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))


def d2(S, K, T, r, sigma):
    return (np.log(S / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))


def delta(S, K, T, r, sigma, option):
    return sy.diff(bachelier(r,K,S,T,sigma,option),S)


def theta(S, K, T, r, sigma, option):
    if option == 'call':
        theta = (-sigma * S * si.norm.pdf(d1(S, K, T, r, sigma))) \
                / (2 * np.sqrt(T)) - r * K * np.exp(-r * T) \
                * si.norm.cdf(d2(S, K, T, r, sigma), 0.0, 1.0)
    if option == 'put':
        theta = (-sigma * S * si.norm.pdf(d1(S, K, T, r, sigma))) \
                / (2 * np.sqrt(T)) + r * K * np.exp(-r * T) \
                * si.norm.cdf(-d2(S, K, T, r, sigma), 0.0, 1.0)

    return theta


def gamma(S, K, T, r, sigma):
    return sy.diff(delta(S, K, T, r, sigma),S)


def vega(S, K, T, r, sigma):
    return S * si.norm.pdf(d1(S, K, T, r, sigma)) * np.sqrt(T)


def rho(S, K, T, r, sigma, option):

    if option == 'call':
        rho = T * K * np.exp(-r * T) * si.norm.cdf(d2(S, K, T, r, sigma), 0.0, 1.0)
    if option == 'put':
        rho = -T * K * np.exp(-r * T) * si.norm.cdf(-d2(S, K, T, r, sigma), 0.0, 1.0)

    return rho

def vanna(S, K, T, r, sigma):
    return np.sqrt(T) * si.norm.pdf(d1(S, K, T, r, sigma)) * (1 - d1(S, K, T, r, sigma))


def volga(S, K, T, r, sigma):
    return np.sqrt(T) * si.norm.pdf(d1(S, K, T, r, sigma)) * (d1(S, K, T, r, sigma)*d2(S, K, T, r, sigma))/sigma

