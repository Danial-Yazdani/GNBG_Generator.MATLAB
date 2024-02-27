import numpy as np
import matplotlib.pyplot as plt


class GNBG:
    def __init__(self):
        np.random.seed(1234)
        self.MaxEvals = 100000
        self.AcceptanceThreshold = 1e-08
        self.Dimension = 5
        self.CompNum = 3
        self.MinCoordinate = -100
        self.MaxCoordinate = 100
        self.CompMinPos = self.initialize_component_position()
        self.CompSigma = self.initialize_component_sigma()
        self.CompH = self.initialize_component_H()

        self.Mu = None
        self.Omega = None
        self.Lambda = None
        self.RotationMatrix = None
        self.OptimumValue = None
        self.OptimumPosition = None
        self.FEhistory = []
        self.FE = 0
        self.AcceptanceReachPoint = np.inf
        self.BestFoundResult = np.inf

    # Initializing the minimum/center position of components
    def initialize_component_position(self): 
        MinRandOptimaPos   = -80
        MaxRandOptimaPos   = 80
        MinExclusiveRange  = -30; # Must be LARGER than GNBG.MinRandOptimaPos
        MaxExclusiveRange  = 30;  # Must be SMALLER than GNBG.MaxRandOptimaPos
        ComponentPositioningMethod = 1 #(1) Random positions with uniform distribution inside the search range
                                       #(2) Random positions with uniform distribution inside a specified range [GNBG.MinRandOptimaPos,GNBG.MaxRandOptimaPos]
                                       #(3) Random positions inside a specified range [GNBG.MinRandOptimaPos,GNBG.MaxRandOptimaPos] but not within the sub-range [GNBG.MinExclusiveRange,GNBG.MaxExclusiveRange]
                                       #(4) Random OVERLAPPING positions with uniform distribution inside a specified range [GNBG.MinRandOptimaPos,GNBG.MaxRandOptimaPos]. Remember to also set GNBG.SigmaPattern to 2. 
        if ComponentPositioningMethod == 1:
           CompMinPos = self.MinCoordinate + (self.MaxCoordinate - self.MinCoordinate) * np.random.rand(self.CompNum, self.Dimension)
        elif ComponentPositioningMethod == 2:
            CompMinPos = MinRandOptimaPos + (MaxRandOptimaPos - MinRandOptimaPos) * np.random.rand(self.CompNum, self.Dimension)
        elif ComponentPositioningMethod == 3:
            lower_range = MinRandOptimaPos + (MinExclusiveRange - MinRandOptimaPos) * np.random.rand(self.CompNum, self.Dimension)  # Generate random numbers in [MinRandOptimaPos, MinExclusiveRange)
            upper_range = MaxExclusiveRange + (MaxRandOptimaPos - MaxExclusiveRange) * np.random.rand(self.CompNum, self.Dimension)  # Generate random numbers in (MaxExclusiveRange, MaxRandOptimaPos]
            selector = np.random.randint(0, 2, size=(self.CompNum, self.Dimension))  # Randomly choose whether to take from lower_range or upper_range
            CompMinPos = (selector * lower_range) + ((1 - selector) * upper_range)
        elif ComponentPositioningMethod == 4:
            CompMinPos = MinRandOptimaPos + np.tile(((MaxRandOptimaPos - MinRandOptimaPos) * np.random.rand(1, self.Dimension)), (self.CompNum, 1))  # Generating o overlapping minimum positions
        else:
            raise ValueError('Warning: Wrong number is chosen for GNBG.ComponentPositioningMethod.')
        return CompMinPos
    

    # Initialize the minimum values of the components
    def initialize_component_sigma(self):
        MinSigma = -99
        MaxSigma = -98
        SigmaPattern = 1 # (1) A random sigma value for EACH component.
                         # (2) A random sigma value for ALL components. It must be used for generating overlapping scenarios, or when the user plans to generate problem instances with multiple global optima.
                         # (3) Manually defined values for sigma.
        if SigmaPattern == 1:
            ComponentSigma = MinSigma + (MaxSigma - MinSigma) * np.random.rand(self.CompNum, 1)
        elif SigmaPattern == 2:
            random_value = MinSigma + (MaxSigma - MinSigma) * np.random.rand()
            ComponentSigma = np.full((self.CompNum, 1), random_value)
        elif SigmaPattern == 3:
            # USER-DEFINED ==> Adjust the size of this array to match the number of components (self.o)
            ComponentSigma = np.array([[-1000], [-950]])
        else:
            raise ValueError('Wrong number is chosen for SigmaPattern.')
        return ComponentSigma
    
    # Defining the elements of diagonal elements of H for components (Configuring condition number)
    def initialize_component_H(self):
        H_pattern = 3  # (1) Condition number is 1 and all elements of principal diagonal of H are set to a user defined value H_value
                       # (2) Condition number is 1 for all components but the elements of principal diagonal of H are different from a component to another and are randomly generated with uniform distribution within the range [Lb_H, Ub_H].
                       # (3) Condition number is random for all components the values of principal diagonal of the matrix H for each component are generated randomly within the range [Lb_H, Ub_H] using a uniform distribution.
                       # (4) Condition number is Ub_H/Lb_H for all components where two randomly selected elements on the principal diagonal of the matrix H are explicitly set to Lb_H and Ub_H. The remaining diagonal elements are generated randomly within the range [Lb_H, Ub_H]. These values follow a Beta distribution characterized by user-defined parameters alpha and beta, where 0 < alpha = beta <= 1.
                       # (5) Condition number is Ub_H/Lb_H for all components where a vector with Dimension equally spaced values between Lb_H and Ub_H is generated. The linspace function is used to create a linearly spaced vector that includes both the minimum and maximum values. For each component, a randomly permutation of this vector is used. 
        Lb_H = 1  # Lower bound for H
        Ub_H = 10**5  # Upper bound for H
        alpha = 0.2  # Example, for Beta distribution
        beta = alpha  # Assuming symmetric distribution        
        if H_pattern == 1:
            H_value = 1
            CompH = H_value * np.ones((self.CompNum, self.Dimension))
        elif H_pattern == 2:
            CompH = (Lb_H + ((Ub_H - Lb_H) * np.random.rand(self.CompNum, 1))) * np.ones((self.CompNum, self.Dimension))
        elif H_pattern == 3:
            CompH = Lb_H + ((Ub_H - Lb_H) * np.random.rand(self.CompNum, self.Dimension))
        elif H_pattern == 4:
            CompH = Lb_H + ((Ub_H - Lb_H) * np.random.beta(alpha, beta, (self.CompNum, self.Dimension)))
            for ii in range(self.CompNum):
                random_indices = np.random.choice(self.Dimension, 2, replace=False)
                CompH[ii, random_indices[0]] = Lb_H
                CompH[ii, random_indices[1]] = Ub_H
        elif H_pattern == 5:
            H_Values = np.linspace(Lb_H, Ub_H, self.Dimension)
            CompH = np.array([np.random.permutation(H_Values) for _ in range(self.CompNum)])
        else:
            raise ValueError('Wrong number is chosen for H_pattern.')
        return CompH
    

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





gnbg = GNBG()


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