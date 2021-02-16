import numpy as np
strike_spreads = []

# set starting parameters
starting_guess = np.array([0.001, 0.5, 0, 0.001])
alpha = len(F) * [starting_guess[0]]
beta = len(F) * [starting_guess[1]]
rho = len(F) * [starting_guess[2]]
nu = len(F) * [starting_guess[3]]
