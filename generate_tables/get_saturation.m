function result = get_saturation( luminance )
coeffs = [0.463467443372641,-0.00353235784321854,0.115988356407610,0.800611655810124,1.73046565180476];
result = ((luminance)*coeffs(1))./sqrt(coeffs(2)+((luminance+coeffs(3))*coeffs(4)).^2)*coeffs(5);