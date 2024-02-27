function chi_pdf_gamma(x, c, gm)
    # gm[3] := gamma(c/2.0)
    # (c -> degree of freedom)
    
    # chi probability density function
    # Arguments:
    # x - value
    # c - degree
    
    pdf = 0.0
    
    if (0 < x < 100)
        y = 0.5 * x * x
        p = 0.5 * c
        
        pdf = exp(-y) * x^(c - 1.0) / (2.0^(p - 1.0) * gm[3])
    end
    
    return pdf
end
