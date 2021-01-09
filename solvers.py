from numba.core.types.containers import Tuple
from numba.core.types.functions import Function
import numpy as np
import json

from numba import jitclass, vectorize, float64

from example_funcs import derivative_lorenz

# Plot in 3D
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

import plotly.express as px
import plotly.graph_objects as go

def read_solver(solver):
	with open("solvers.json", "r") as file:
		solvers = json.load(file)

	BT = ButcherTableau(**{att: np.array(eval(expression)) for att, expression in solvers['explicit_methods'][solver].items()}, name = solver)
	return BT

class ButcherTableau(object):
	def __init__(self, alpha: np.array, beta: np.array, gamma: np.array, name = "Custom Solver"):
		self.alpha = alpha
		self.beta = beta
		self.gamma = gamma

		assert((len(alpha) == len(beta)))
		assert(gamma.shape[0] == gamma.shape[1])

	def __str__(self):
		rep = ""
		if self.title is not None: rep += "Model: " + self.title + "\n\n"

		# rep += "Butcher Tableau:\n===============\n\n"
		# for b in self.beta:
		# 	rep += str(b) + '|' + [str()]

	def __repr__(self):
		return self.title

# @jitclass([
# 	('alpha', float64[:]),
# 	('beta', float64[:]),
# 	('gamma', float64[:,:]),
# 	('positions', float64[:,:]),
# 	# ('behaviour', Function),
# 	# ('extra_args', Tuple),
# 	('clock', float),
# 	('ttl', float64[:])
# ])
class Swarm(object):
	def __init__(self, tableau: ButcherTableau, function: Function, *args, size: int = 1000, dimensions = 2, initial_t = 0, ):
		self.alpha, self.beta, self.gamma = tableau.alpha, tableau.beta, tableau.gamma
		self.tableau = tableau

		np.random.seed(1)

		self.positions = np.random.rand(dimensions, size)*10

		# self.particles = [Particle(self.positions[:,i], ttl = 100) for i in range(size)]

		self.behaviour = function
		self.extra_args = args

		self.clock = initial_t

		self.ttl = np.empty(size)

	def swarm_step(self, step: float = 0.1) -> float:

		# Set initial stepsize
		h = step

		# Preallocate array
		derivative_evaluations = np.zeros([self.positions.shape[0], self.positions.shape[1], len(self.beta)])

		# Calculate derivative 
		for j in range(len(self.beta)): # j is index of particular derivative
			derivative_evaluations[:,:,j] = self.behaviour(self.clock + self.beta[j]*h, self.positions + h*np.matmul(derivative_evaluations, self.gamma[j,:]), *self.extra_args)

		# Calculate gradient at current tn for each embedded method
		f = np.matmul(derivative_evaluations, self.alpha.T)

		# Calculate estimates for y for each method in embedded method
		y_estimates =  h*f + self.positions

		# Get estimate of total error
		largest_error = np.max(abs(np.diff(y_estimates, axis = 0)))
				
		self.clock += h
		self.positions = y_estimates

class Particle(Swarm):
	def __init__(self, position: np.ndarray, ttl: float) -> None:
		self.position = position
		self.ttl = ttl
		
	def calculate_step():
		pass


if __name__ == "__main__":
	BT = read_solver("RK4")
	sigma = 10.
	rho = 28.
	beta = 8. / 3.
	args = sigma, rho, beta
	particles = Swarm(BT, derivative_lorenz, size = 50, dimensions = 3, *args)

	# def animated(d):
	# 	particles.swarm_step()
	# 	ax.scatter(particles.positions[0,:], particles.positions[1,:], particles.positions[2,:])

	# animated_step = animated
	# fig = plt.figure()
	# fig.subplots_adjust(left=0, bottom=0, right=1, top=1)
	# ax = fig.add_subplot(111, projection='3d')
	# ax.set_facecolor((0.5, 0.5, 0.5))

	# animation = FuncAnimation(fig, animated, interval = 300)
	# plt.show()

	fig = px.scatter_3d()
	
	for _ in range(10):
		particles.swarm_step()
		fig.add_trace(
			go.Scatter3d(
				x=particles.positions[0,:],
				y=particles.positions[1,:],
				z=particles.positions[2,:],
				mode = 'markers',
				showlegend=False)
		)

	fig.show()

	