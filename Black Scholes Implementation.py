# Black Scholes Implementation (Price of Call Option)

import numpy as np
from scipy.stats import norm

# Parameters

r = 0.05
K = 1
t = 2
T = 5
sigma = np.sqrt(2*r)
x = 5

#Black Scholes Formula

def Black_Scholes(r, K, T, sigma, x):
    d1 = (np.log(x/K) + (r + sigma**2/2)*(T-t))/(sigma*np.sqrt(T-t))
    d2 = d1 - sigma*np.sqrt(T-t)
    price = x*norm.cdf(d1) - K*np.exp(-r*(T-t))*norm.cdf(d2)
    return price

print("Price of call option is:", round(Black_Scholes(r, K, T, sigma, x), 2))