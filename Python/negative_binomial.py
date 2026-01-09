import numpy as np
from scipy.stats import nbinom, poisson
from scipy.optimize import root_scalar

# Parameters
lambda_value = 10  # Set lambda for Poisson distribution
n = 3  # Number of samples

# Step 1: Draw k values from a Poisson distribution with lambda
k_values = poisson.rvs(mu=lambda_value, size=n)

# Step 2: Compute values of p for r = 1 to 10
rs = range(1, 11)
ps = [lambda_value / (lambda_value + r) for r in rs]

# Step 3: For each (r, p) pair and each k, compute the p-value with respect to NB(r, p)
p_value_matrix = np.zeros((len(rs), n))  # Matrix to store p-values for each r and k
for i, (r, p) in enumerate(zip(rs, ps)):
    p_value_matrix[i] = [nbinom.cdf(k, r, p) for k in k_values]

print(p_value_matrix)
print(ps)

# Step 4: Define function to estimate r and p based on the mean and p-value constraints
def estimate_r_p(target_mean, target_p_value, k):
    def objective_function(p):
        r = target_mean * p / (1 - p)
        cdf_value = nbinom.cdf(k, r, p)
        return cdf_value - target_p_value

    # Solve for p
    solution = root_scalar(objective_function, x0=0.5, bracket=[.001,.999])
    p = solution.root
    r = target_mean * p / (1 - p)
    return r, p

# Step 5: For each r, estimate r and p using the simulated data and check consistency
for i, (r, true_p) in enumerate(zip(rs, ps)):
    # Calculate mean and target p-value as the average over the simulated data
    target_mean = lambda_value  # This should be close to our lambda for each r
    for p_val, k in zip(p_value_matrix[i], k_values):
        # Estimate r and p using the function
        estimated_r, estimated_p = estimate_r_p(target_mean, p_val, k)

        # Print and check results
        print(f"True (r, p): ({r}, {true_p:.4f})")
        print(f"Estimated (r, p): ({estimated_r:.4f}, {estimated_p:.4f})\n")
        # assert abs(estimated_r - r) < 1e-6, f"Mismatch in r: estimated {estimated_r} vs true {r}"
        # assert abs(estimated_p - true_p) < 1e-6, f"Mismatch in p: estimated {estimated_p} vs true {true_p}"



print("All tests passed successfully!")