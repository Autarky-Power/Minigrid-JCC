# ======================================
# PROVIDE Data and Sampling
# ======================================

# Set dimension of each probability constraint
dim_x = length(T) + outage
dim_cons = outage + 1
SampleSize = 8000

# Set probability level
p = inputs.p

# setup mean
μ = inputs.μ

# setup covariance
Σ = inputs.Σ

# Compute cholesky of covariance
L = cholesky(Σ, check=true).L

# Compute Sample on the sphere
SampleOnSphere = randn(SampleSize, dim_cons*(dim_x-dim_cons+1))

for j in 1:dim_x-dim_cons+1
    i = dim_cons * (j-1) + 1
    for n = 1:SampleSize
    # normalize sample (projection to sphere)
        v = SampleOnSphere[n, i:i+dim_cons-1]'
        v = v / norm(v)

        # transformation with cholesky of covariance
        v = L[j:j+dim_cons-1, j:j+dim_cons-1] * v'

        # save result
        SampleOnSphere[n, i:i+dim_cons-1] = v'
    end 
end
