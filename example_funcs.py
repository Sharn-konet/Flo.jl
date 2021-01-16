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

		fx = alpha*(y-x)
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

	def shimizu_morioka(t,v, alpha: float = 0.75, beta: float = 0.45):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = y
		fy = (1-z)*x - alpha*y
		fz = x*x - beta*z

		return np.array([fx, fy, fz])

	def sakarya(t, v, alpha: float = 0.4, beta: float = 0.3):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = -x + y + y*z
		fy = -x-y + alpha*x*z
		fz = z - beta*x*y

		return np.array([fx, fy, fz])

	def rucklidge(t, v, kappa: float = 2, alpha: float = 6.7):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = -kappa*x + alpha*y - y*z
		fy = x
		fz = -z + y*y

		return np.array([fx, fy, fz])

	def rossler(t, v, alpha: float = 0.2, beta: float = 0.2, sigma: float = 5.7):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = -(y+z)
		fy = x + alpha*y
		fz = beta + z*(x - sigma)

		return np.array([fx, fy, fz])

	def rayleigh_benard(t, v, alpha: float = 9, gamma: float = 12, beta: float = 5):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = -alpha*x + alpha*y
		fy = gamma*x - y - x*z
		fz = x*y - beta*z

		return np.array([fx, fy, fz])

	def qi_chen(t, v, alpha: float = 38, beta: float = 8/3, sigma: float = 80):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = alpha * (y-x) +  y*z
		fy = sigma*x + y - x*z
		fz = x*y - beta*z

		return np.array([fx, fy, fz])

	def qi(t, v, alpha: float = 30, beta: float = 10, sigma: float = 1, delta: float = 10):

		x, y, z, w = v[0,:], v[1,:], v[2,:], v[3,:]

		fx = alpha*(y-x) + y*z*w
		fy = beta*(x+y) - x*z*w
		fz = -sigma*z + x*y*w
		fw = -delta*w + x*y*z

		return np.array([fx, fy, fz, fw])

	def nose_hoover(t, v, alpha: float = 1.5):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = y
		fy = -x + y*z
		fz = alpha-y*y

		return np.array([fx, fy, fz])

	def newton_leipnik(t, v, alpha: float = 0.4, beta: float = 0.175):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = -alpha*x + y+ 10*y*z
		fy = -x-0.4*y + 5*x*z
		fz = beta*z-5*x*y

		return np.array([fx, fy, fz])

	def lü_chen(t, v, alpha: float = -10, beta: float = -4, sigma: float = 18.1):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = -((alpha*beta)/(alpha+beta))*x -y*z + sigma
		fy = alpha*y + x*z
		fz = beta*z + x*y

		return np.array([fx, fy, fz])

	def aizawa(t, v, alpha: float = 0.95, beta: float = 0.7, gamma:float = 0.6, delta: float = 3.5, epsilon: float = 0.25, zeta: float = 0.1):

		x, y, z = v[0,:], v[1,:], v[2,:]

		fx = (z-beta)*x-y
		fy = delta*x+(z-beta)*y
		fz = gamma + alpha*z -(1/3)*z*z*z - (x*x + y*y)*(1 + epsilon*z) + zeta*z*x*x*x

		return np.array([fx, fy, fz])