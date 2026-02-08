import numpy as np

num_samples = 3000
mean = 1.0
three_sigma_percentage = 0.1 # 10% variation
sigma = (mean * three_sigma_percentage)/3

seeds = np.random.normal(mean, sigma, num_samples)

filename = 'accel_seeds.dat'
with open(filename, 'w') as f:
    for s in seeds:
        formatted_val = f"{s:15.7e}"
        f.write(f"{formatted_val}\n")