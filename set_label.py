from set_swaption import *
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