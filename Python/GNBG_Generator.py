import os
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt


class GNBG:
    def __init__(self, MaxEvals, AcceptanceThreshold, Dimension, CompNum, MinCoordinate, MaxCoordinate, CompMinPos, CompSigma, CompH, Mu, Omega, Lambda, RotationMatrix, OptimumValue, OptimumPosition):
        self.MaxEvals = MaxEvals
        self.AcceptanceThreshold = AcceptanceThreshold
        self.Dimension = Dimension
        self.CompNum = CompNum
        self.MinCoordinate = MinCoordinate
        self.MaxCoordinate = MaxCoordinate
        self.CompMinPos = CompMinPos
        self.CompSigma = CompSigma
        self.CompH = CompH
        self.Mu = Mu
        self.Omega = Omega
        self.Lambda = Lambda
        self.RotationMatrix = RotationMatrix
        self.OptimumValue = OptimumValue
        self.OptimumPosition = OptimumPosition
        self.FEhistory = []
        self.FE = 0
        self.AcceptanceReachPoint = np.inf
        self.BestFoundResult = np.inf

    def fitness(self, X):
        SolutionNumber = X.shape[0]
        result = np.nan * np.ones(SolutionNumber)
        for jj in range(SolutionNumber):
            x = X[jj, :].reshape(-1, 1)  # Ensure column vector
            f = np.nan * np.ones(self.CompNum)
            for k in range(self.CompNum):
                if len(self.RotationMatrix.shape) == 3:
                    rotation_matrix = self.RotationMatrix[:, :, k]
                else:
                    rotation_matrix = self.RotationMatrix

                a = self.transform((x - self.CompMinPos[k, :].reshape(-1, 1)).T @ rotation_matrix.T, self.Mu[k, :], self.Omega[k, :])
                b = self.transform(rotation_matrix @ (x - self.CompMinPos[k, :].reshape(-1, 1)), self.Mu[k, :], self.Omega[k, :])
                f[k] = self.CompSigma[k] + (a @ np.diag(self.CompH[k, :]) @ b) ** self.Lambda[k]

            result[jj] = np.min(f)
            if self.FE > (self.MaxEvals-1):
                return result
            self.FE += 1
            self.FEhistory = np.append(self.FEhistory, result[jj])
            if self.BestFoundResult > result[jj]:
                self.BestFoundResult = result[jj]
            if abs(self.FEhistory[self.FE-1] - self.OptimumValue) < self.AcceptanceThreshold and np.isinf(self.AcceptanceReachPoint):
                self.AcceptanceReachPoint = self.FE
        return result

    def transform(self, X, Alpha, Beta):
        Y = X.copy()
        tmp = (X > 0)
        Y[tmp] = np.log(X[tmp])
        Y[tmp] = np.exp(Y[tmp] + Alpha[0] * (np.sin(Beta[0] * Y[tmp]) + np.sin(Beta[1] * Y[tmp])))
        tmp = (X < 0)
        Y[tmp] = np.log(-X[tmp])
        Y[tmp] = -np.exp(Y[tmp] + Alpha[1] * (np.sin(Beta[2] * Y[tmp]) + np.sin(Beta[3] * Y[tmp])))
        return Y


# Get the current script's directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Define the path to the folder where you want to read/write files
folder_path = os.path.join(current_dir)

# Initialization
ProblemIndex = 1  # Choose a problem instance range from f1 to f24

# Check if the chosen problem index is within the valid range
if ProblemIndex < 1 or ProblemIndex > 24:
    raise ValueError("Invalid problem index. Please choose a number from 1 to 24.")

# Preparation and loading of the GNBG parameters based on the chosen problem instance
if 1 <= ProblemIndex <= 24:
    filename = f'f{ProblemIndex}.mat'
    GNBG_tmp = loadmat(os.path.join(folder_path, filename))['GNBG']
    MaxEvals = np.array([item[0] for item in GNBG_tmp['MaxEvals'].flatten()])[0, 0]
    AcceptanceThreshold = np.array([item[0] for item in GNBG_tmp['AcceptanceThreshold'].flatten()])[0, 0]
    Dimension = np.array([item[0] for item in GNBG_tmp['Dimension'].flatten()])[0, 0]
    CompNum = np.array([item[0] for item in GNBG_tmp['o'].flatten()])[0, 0]  # Number of components
    MinCoordinate = np.array([item[0] for item in GNBG_tmp['MinCoordinate'].flatten()])[0, 0]
    MaxCoordinate = np.array([item[0] for item in GNBG_tmp['MaxCoordinate'].flatten()])[0, 0]
    CompMinPos = np.array(GNBG_tmp['Component_MinimumPosition'][0, 0])
    CompSigma = np.array(GNBG_tmp['ComponentSigma'][0, 0], dtype=np.float64)
    CompH = np.array(GNBG_tmp['Component_H'][0, 0])
    Mu = np.array(GNBG_tmp['Mu'][0, 0])
    Omega = np.array(GNBG_tmp['Omega'][0, 0])
    Lambda = np.array(GNBG_tmp['lambda'][0, 0])
    RotationMatrix = np.array(GNBG_tmp['RotationMatrix'][0, 0])
    OptimumValue = np.array([item[0] for item in GNBG_tmp['OptimumValue'].flatten()])[0, 0]
    OptimumPosition = np.array(GNBG_tmp['OptimumPosition'][0, 0])
else:
    raise ValueError('ProblemIndex must be between 1 and 24.')

gnbg = GNBG(MaxEvals, AcceptanceThreshold, Dimension, CompNum, MinCoordinate, MaxCoordinate, CompMinPos, CompSigma, CompH, Mu, Omega, Lambda, RotationMatrix, OptimumValue, OptimumPosition)

# Set a random seed for the optimizer
np.random.seed()  # This uses a system-based source to seed the random number generator


#Your optimization algorithm goes here
X = np.random.rand(10000, Dimension) # This is for generating a random population of 5 solutions for testing GNBG
# Calculating the fitness=objective values of the population using the GNBG function. The result is a 5x1 vector of objective values
result = gnbg.fitness(X) 

# After running the algorithm, the best fitness value is stored in gnbg.BestFoundResult
# The function evaluation number where the algorithm reached the acceptance threshold is stored in gnbg.AcceptanceReachPoint
# For visualizing the convergence behavior, the history of the objective values is stored in gnbg.FEhistory, however it needs to be processed as follows:

convergence = []
best_error = float('inf')
for value in gnbg.FEhistory:
    error = abs(value - OptimumValue)
    if error < best_error:
        best_error = error
    convergence.append(best_error)

# Plotting the convergence
plt.plot(range(1, len(convergence) + 1), convergence)
plt.xlabel('Function Evaluation Number (FE)')
plt.ylabel('Error')
plt.title('Convergence Plot')
#plt.yscale('log')  # Set y-axis to logarithmic scale
plt.show()