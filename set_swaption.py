from input_output import*
import numpy as np
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