import numpy as np

def derivative_lorenz(t, y, sigma, rho, beta):
	"""
	Compute the derivatives for the Lorenz system.

	Args:
		t (float): independent variable i.e. time (s).
		y (ndarray): current x, y and z in Lorenz system.
		sigma (float): system parameter of Lorenz system.
		rho (float): system parameter of Lorenz system.
		beta (float): system parameter of Lorenz system.

	Returns:
		f (ndarray): derivatives of x, y and z in Lorenz system.
	"""
	
	fx = sigma*(y[1,:] - y[0,:])
	fy = y[0,:] * (rho - y[2,:]) - y[1,:]
	fz = y[0,:] * y[1,:] - beta*y[2,:]

	return np.array([fx, fy, fz])