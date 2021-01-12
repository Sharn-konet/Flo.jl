import numpy as np

class StrangeAttractors:

	def lorenz(t, v, sigma: float = 10, rho: float = 28, beta: float = 8/3):
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

		x, y, z = v[0,:], v[1,:], v[2,:]
		
		fx = sigma*(y - x)
		fy = x * (rho - z) - y
		fz = x * y - beta*z

		return np.array([fx, fy, fz])

	def TSUCS1(t, v, alpha: float = 40, beta: float = 0.5 + 1/3, delta: float = 0.5, epsilon: float = 0.65, zeta: float = 20):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = alpha*(y-x) + delta*x*z
		fy = zeta*y -x*z
		fz = beta*z + x*y - epsilon*x*x

		return np.array([fx, fy, fz])


	def TSUCS2(t, v, a: float = 40., b: float = 55, c: float = 1.5 + 1/3, d: float = 0.16, e: float = 0.65, f: float = 20.0):
		
		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = a*(y - x) + d*x*z
		fy = (b*x) - (x*z) + f*y
		fz = (c*z) + (x * y) - (e + x*x)

		return np.array([fx, fy, fz])

	def yuwang(t, v, alpha: float = 10, beta: float = 40, sigma: float = 2, delta: float = 2.5):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = alpha(y-x)
		fy = beta*x - sigma*x*z
		fz = np.exp(x*y) - delta*z

		return np.array([fx, fy, fz])

	def wimol_banlue(t, v, alpha: float = 2):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = y-x
		fy = -z*np.tanh(x)
		fz = -alpha * x*y + np.absolute(y)

		return np.array([fx, fy, fz])
		
	def wang_sun(t, v, alpha: float = 0.2, beta:float = -0.01, sigma: float = 1.0, delta: float = -0.4, epsilon: float = -1., zeta: float = -1.):
		
		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = alpha*x + sigma*y*z
		fy = beta*x + delta*y -x*z
		fz = epsilon*z + zeta*x*y

		return np.array([fx, fy, fz])

	def thomas(t, v, beta: float = 0.19):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = beta*x + np.sin(y)
		fy = -beta*y + np.sin(z)
		fz = -beta*z + np.sin(x)

		return np.array([fx, fy, fz])
