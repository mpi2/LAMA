from numpy import sqrt, log

def sigma2Gamma(sigma):
    '''Function to convert standard deviation (sigma) to FWHM (Gamma)'''
    return sigma * sqrt(2 * log(2)) * 2 / sqrt(2)

Gamma = sigma2Gamma(10)

# prints 16.651092223153956, which is what you saw in your graph

def Gamma2sigma(Gamma):
    '''Function to convert FWHM (Gamma) to standard deviation (sigma)'''
    return Gamma * sqrt(2) / ( sqrt(2 * log(2)) * 2 )

fwhm = 200.0
pixel_size = 28.0
pixel_width = fwhm /pixel_size
sigma = Gamma2sigma(pixel_width)
print sigma