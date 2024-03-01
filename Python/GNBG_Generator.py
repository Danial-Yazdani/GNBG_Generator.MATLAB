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
        self.Mu, self.Omega = self.initialize_modality_parameters()
        self.Lambda = self.initialize_lambda()
        self.RotationMatrix = self.define_rotation_matrices()
        self.OptimumValue = np.min(self.CompSigma)
        self.OptimumPosition = self.CompMinPos[np.argmin(self.CompSigma), :]
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
            raise ValueError('Warning: Wrong number is chosen for ComponentPositioningMethod.')
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
    
    # Defining the parameters used in transformation function (Configuring modality of components and their basin local optima)
    def initialize_modality_parameters(self):
        MinMu = 0.2
        MaxMu = 0.5
        MinOmega = 5
        MaxOmega = 50
        localModalitySymmetry = 3  # (1) Unimodal, smooth, and regular components
                                   # (2) Multimodal symmetric components
                                   # (3) Multimodal asymmetric components
                                   # (4) Manually defined values
        if localModalitySymmetry == 1:
            Mu = np.zeros((self.CompNum, 2))
            Omega = np.zeros((self.CompNum, 4))
        elif localModalitySymmetry == 2:
            Mu = np.tile(MinMu + (MaxMu - MinMu) * np.random.rand(self.CompNum, 1), (1, 2))
            Omega = np.tile(MinOmega + (MaxOmega - MinOmega) * np.random.rand(self.CompNum, 1), (1, 4))
        elif localModalitySymmetry == 3:
            Mu = MinMu + (MaxMu - MinMu) * np.random.rand(self.CompNum, 2)
            Omega = MinOmega + (MaxOmega - MinOmega) * np.random.rand(self.CompNum, 4)
        elif localModalitySymmetry == 4:
            # Assuming self.CompNum is defined and matches the required shape
            # User-defined values; adjust sizes as needed
            Mu = np.array([[1, 1]] * self.CompNum)  # Adjust size as necessary
            Omega = np.array([[10, 10, 10, 10]] * self.CompNum)  # Adjust size as necessary
        else:
            raise ValueError('Wrong number is chosen for localModalitySymmetry.')
        return Mu, Omega
    
    # Defining the linearity of the basin of components
    def initialize_lambda(self):
        MaxLambda = 1
        MinLambda = 1
        LambdaValue4ALL = 0.25  # Default value; adjust as needed
        LambdaConfigMethod = 1  # 1 All lambda are set to LambdaValue4ALL
                                # 2 Randomly set lambda of each component in [MinLambda,MaxLambda]. Note that large ranges may result in existence of invisible components
        if LambdaConfigMethod == 1:
            Lambda = np.full((self.CompNum, 1), LambdaValue4ALL)
        elif LambdaConfigMethod == 2:
            Lambda = MinLambda + (MaxLambda - MinLambda) * np.random.rand(self.CompNum, 1)
        else:
            raise ValueError('Wrong number is chosen for LambdaConfigMethod.')
        return Lambda
    
    def define_rotation_matrices(self):
        MinAngle = -np.pi
        MaxAngle = np.pi
        Rotation = 2  # (1) Without rotation
                      # (2) Random Rotation for all components==> Fully Connected Interaction with random angles for each plane of each component
                      # (3) Random Rotation for all components==>For each component, interactions are defined randomly based on connection probability threshold GNBG.ConnectionProbability and random angles. GNBG.ConnectionProbability can be different for each component.
                      # (4) Rotation based on a random Angle for each component (fully connected with an angle for all planes in each component)
                      # (5) Rotation based on the specified Angle for all components (fully connected with an angle for all planes in all component)
                      # (6) Rotation based on the random Angles to generate chain-like variable interaction structure
                      # (7) Generating partially separable variable interaction structure with user defined sizes and angles for each group of variables
        RotationMatrix = np.nan * np.zeros((self.CompNum, self.Dimension, self.Dimension))
        if Rotation == 1:
            RotationMatrix = np.array([np.eye(self.Dimension) for _ in range(self.CompNum)])
        elif Rotation in [2, 3, 4, 5, 6, 7]:
            if Rotation == 2:
                for ii in range(self.CompNum):
                    ThetaMatrix = np.zeros((self.Dimension, self.Dimension))
                    upper_triangle_indices = np.triu_indices(self.Dimension, k=1)
                    ThetaMatrix[upper_triangle_indices] = MinAngle + (MaxAngle - MinAngle) * np.random.rand(len(upper_triangle_indices[0]))
                    RotationMatrix[ii, :, :] = self.rotation(ThetaMatrix)
            elif Rotation == 3:
                MinConProb = 0.5
                MaxConProb = 0.75
                ConnectionProbability = MinConProb + (MaxConProb - MinConProb) * np.random.rand(self.CompNum, 1)
                for ii in range(self.CompNum):
                    ThetaMatrix = np.zeros((self.Dimension, self.Dimension))
                    for i in range(self.Dimension):
                        for j in range(i + 1, self.Dimension):
                            if np.random.rand() < ConnectionProbability[ii]:
                                ThetaMatrix[i, j] = MinAngle + (MaxAngle - MinAngle) * np.random.rand()
                    RotationMatrix[ii, :, :] = self.rotation(ThetaMatrix)
            elif Rotation == 4:
                randomAngle = MinAngle + (MaxAngle - MinAngle) * np.random.rand(self.CompNum, 1)
                for ii in range(self.CompNum):
                    ThetaMatrix = np.full((self.Dimension, self.Dimension), randomAngle[ii])
                    lower_triangle_indices = np.tril_indices(ThetaMatrix.shape[0])
                    ThetaMatrix[lower_triangle_indices] = 0 # Set all elements on and below the diagonal to 0
                    RotationMatrix[ii, :, :] = self.rotation(ThetaMatrix)
            elif Rotation == 5:
                for ii in range(self.CompNum):
                    SpecificAngles = 1.48353
                    ThetaMatrix = np.full((self.Dimension, self.Dimension), SpecificAngles)
                    lower_triangle_indices = np.tril_indices(ThetaMatrix.shape[0])
                    ThetaMatrix[lower_triangle_indices] = 0 # Set all elements on and below the diagonal to 0
                    RotationMatrix[ii, :, :] = self.rotation(ThetaMatrix)
            elif Rotation == 6:
                for ii in range(self.CompNum):
                    ThetaMatrix = np.zeros((self.Dimension, self.Dimension))
                    for i in range(self.Dimension - 1):
                        ThetaMatrix[i, i + 1] = MinAngle + (MaxAngle - MinAngle) * np.random.rand()
                    RotationMatrix[ii, :, :] = self.rotation(ThetaMatrix)
            elif Rotation == 7:
                S = np.array([3, 4, 3])  # Number of variables in each group
                Theta = [np.pi/4, 3*np.pi/4, np.pi/8]  # Angle for each group
                if sum(S) > self.Dimension or len(S) != len(Theta):
                    raise ValueError("The sum of elements in S exceeds Dimension, or the sizes of S and Theta are not equal.")
                for ii in range(self.CompNum):
                    ThetaMatrix = np.zeros((self.Dimension, self.Dimension))
                    allVars = np.random.permutation(self.Dimension)
                    groupStart = 0
                    for jj, size in enumerate(S):
                        groupEnd = groupStart + size
                        groupVars = allVars[groupStart:groupEnd]
                        for var1 in groupVars:
                            for var2 in groupVars:
                                if var1 < var2:  # Only for elements above the diagonal
                                    ThetaMatrix[var1, var2] = Theta[jj]
                        groupStart = groupEnd
                    RotationMatrix[ii, :, :] = self.rotation(ThetaMatrix)            
        else:
            raise ValueError('Wrong number is chosen for Rotation.')
        return RotationMatrix

    def rotation(self, teta):
        R = np.eye(self.Dimension)
        for p in range(self.Dimension - 1):
            for q in range(p + 1, self.Dimension):
                if teta[p, q] != 0:
                    G = np.eye(self.Dimension)
                    cos_val = np.cos(teta[p, q])
                    sin_val = np.sin(teta[p, q])
                    G[p, p], G[q, q] = cos_val, cos_val
                    G[p, q], G[q, p] = -sin_val, sin_val
                    R = np.dot(R, G)
        return R

    def fitness(self, X):
        SolutionNumber = X.shape[0]
        result = np.nan * np.ones(SolutionNumber)
        for jj in range(SolutionNumber):
            x = X[jj, :].reshape(-1, 1)  # Ensure column vector
            f = np.nan * np.ones(self.CompNum)
            for k in range(self.CompNum):
                if len(self.RotationMatrix.shape) == 3:
                    rotation_matrix = self.RotationMatrix[k, :, :]
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
X = np.random.rand(10000, gnbg.Dimension) # This is for generating a random population of 5 solutions for testing GNBG
# Calculating the fitness=objective values of the population using the GNBG function. The result is a 5x1 vector of objective values
result = gnbg.fitness(X) 

# After running the algorithm, the best fitness value is stored in gnbg.BestFoundResult
# The function evaluation number where the algorithm reached the acceptance threshold is stored in gnbg.AcceptanceReachPoint
# For visualizing the convergence behavior, the history of the objective values is stored in gnbg.FEhistory, however it needs to be processed as follows:

convergence = []
best_error = float('inf')
for value in gnbg.FEhistory:
    error = abs(value - gnbg.OptimumValue)
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