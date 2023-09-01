"""
    ProbFunction(x, SampleOnSphere, mu)

This function computes the probability of a given vector `x` 
being less than or equal to a normal random variable `xi` 
with mean `mu` and a given covariance matrix, based on a sample
of points distributed on the unit sphere.
"""
function ProbFunction(x, SampleOnSphere, mu)
    
    #   Phi(x) :=  P ( xi <= x )
    #
    #        x  : given vector
    #        xi : normal random variable ~ N(mu,Σ)
    #
    
    # Phi(x) computation
    
    SampleSize, M = size(SampleOnSphere)
    PROB = 0
    
    # Required gamma values for Chi-distibution
    gm = [loggamma(M/2.0), loggamma(M/2.0 + 1.0), gamma(M/2.0)]
    
    
    # setup gradient
    grad = zeros(size(x, 1), 1)
    
    for sample = 1:SampleSize
        
        # Select sample on sphere
        v = SampleOnSphere[sample, :]'
        
        # Compute radius r (Inner problem)
        
        # Pre-Initialization
        r = 1e99
        active = 0
        
        y = x - mu
        
        if maximum(y) < 0
            error("ProbFunction() ... current x no slater point")
        end
        
        for k = 1:M
            
            # k-th component
            rstep = y[k] / v[k]
            
            if rstep > 0
                if rstep < r
                    r = rstep
                    active = k
                end
            end
            
        end
        
        # Compute chi-distribution of r
        
        # PROB
        if r > 0
            chi = chi_cdf_gamma(r, M, gm)
            PROB += chi
        end
        
        # Gradient of ProbFunction
        
       
            
        if r > 0
            
            if active > 0
                
                # active inequality for r
                gradz = zeros(M, 1)
                gradz[active] = 1
                factor = dot(gradz', v)
                factor = chi_pdf_gamma(r, M, gm) / factor
                
                # update grad_vector components
                gradx = zeros(M, 1)
                gradx[active] = -1
                grad -= factor .* gradx
                
            end
            
        end
              
    end
    
    PROB = PROB / SampleSize
    
    # Value of constraint:  probLevel - Phi(u) <= 0
    prob = PROB
    
    # Gradient of constraint: - Grad(Phi(u))
    grad = grad / SampleSize
        
    return (prob, grad)
end