%***********************************GNBG***********************************
%Author: Danial Yazdani
%Last Edited: December 14, 2023
%Title: Generalized Numerical Benchmark Generator (GNBG)
% --------
%Description: 
%          This function is the GNBG problem instance generator. By
%          manupulating the parameters in this code, users can generate a
%          variaty of problem instances.
% --------
%Refrence: 
%           D. Yazdani, M. N. Omidvar, D. Yazdani, K. Deb, and A. H. Gandomi, "GNBG: A Generalized
%           and Configurable Benchmark Generator for Continuous Numerical Optimization," arXiv prepring	arXiv:2312.07083, 2023.
% 
%If you are using GNBG and this code in your work, you should cite the reference provided above.       
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Danial Yazdani
% e-mail: danial DOT yazdani AT gmail DOT com
% Copyright notice: (c) 2023 Danial Yazdani
%************************************************************************** 
function [GNBG] = BenchmarkGenerator
GNBG                     = [];
GNBG.FirstPoint          = 100000;%Sampling the best found solution's objective value when GNBG.FirstPoint function evaluations has been reached.
GNBG.SecondPoint         = 250000;%Sampling the best found solution's objective value when GNBG.SecondPoint function evaluations has been reached.
GNBG.MaxEvals            = 500000;%Maximum function evaluation number
GNBG.AcceptanceThreshold = 1e-08;%The number of function evaluations used to find a solution whose error is GNBG.AcceptanceThreshold, will be saved as a performance indicator.
GNBG.Dimension           = 30;%Dimension number
GNBG.o                   = 1;%Number of components
%% Initializing the minimum/center position of components
GNBG.MinCoordinate       = -100;
GNBG.MaxCoordinate       = 100;
GNBG.MinRandOptimaPos    = -80;
GNBG.MaxRandOptimaPos    = 80;
GNBG.MinExclusiveRange   = -30;%Must be LARGER than GNBG.MinRandOptimaPos
GNBG.MaxExclusiveRange   = 30;%Must be SMALLER than GNBG.MaxRandOptimaPos
GNBG.ComponentPositioningMethod =4;%(1) Random positions with uniform distribution inside the search range
                                   %(2) Random positions with uniform distribution inside a spacefied range [GNBG.MinRandOptimaPos,GNBG.MaxRandOptimaPos]
                                   %(3) Random positions inside a spacefied range [GNBG.MinRandOptimaPos,GNBG.MaxRandOptimaPos] but not within the sub-range [GNBG.MinExclusiveRange,GNBG.MaxExclusiveRange]
                                   %(4) Random OVERLAPPING positions with uniform distribution inside a spacefied range [GNBG.MinRandOptimaPos,GNBG.MaxRandOptimaPos]. Remember to also set GNBG.SigmaPattern to 2. 
switch GNBG.ComponentPositioningMethod
    case 1
        GNBG.Component_MinimumPosition  = GNBG.MinCoordinate + (GNBG.MaxCoordinate-GNBG.MinCoordinate)*rand(GNBG.o,GNBG.Dimension);
    case 2
        GNBG.Component_MinimumPosition  = GNBG.MinRandOptimaPos + (GNBG.MaxRandOptimaPos-GNBG.MinRandOptimaPos)*rand(GNBG.o,GNBG.Dimension);
    case 3
        lower_range = GNBG.MinRandOptimaPos + (GNBG.MinExclusiveRange - GNBG.MinRandOptimaPos) * rand(GNBG.o, GNBG.Dimension);% Generate random numbers in [GNBG.MinRandOptimaPos, GNBG.MinExclusiveRange)
        upper_range = GNBG.MaxExclusiveRange + (GNBG.MaxRandOptimaPos - GNBG.MaxExclusiveRange) * rand(GNBG.o, GNBG.Dimension);% Generate random numbers in (GNBG.MaxExclusiveRange, GNBG.MaxRandOptimaPos]
        selector = randi([0, 1], GNBG.o, GNBG.Dimension);% Randomly choose whether to take from lower_range or upper_range
        GNBG.Component_MinimumPosition = (selector .* lower_range) + ((1 - selector) .* upper_range);
    case 4
        GNBG.Component_MinimumPosition  = GNBG.MinRandOptimaPos + repmat(((GNBG.MaxRandOptimaPos-GNBG.MinRandOptimaPos)*rand(1,GNBG.Dimension)),GNBG.o,1);%Generating GNBG.o overlapping minimum positions
    otherwise
        warning('Wrong number is chosen for GNBG.ComponentPositioningMethod.')
end
%% Defining the minimum values of the components
GNBG.MinSigma       = -99;
GNBG.MaxSigma       = -98;
GNBG.SigmaPattern   = 1;%(1) A random sigma value for EACH component.
                        %(2) A random sigma value for ALL components. It must be used for generating overlapping scenarios, or when the user plans to generate problem instances with multiple global optima.
                        %(3) Manually defined values for sigma.
switch GNBG.SigmaPattern
    case 1
        GNBG.ComponentSigma = GNBG.MinSigma + (GNBG.MaxSigma-GNBG.MinSigma)*rand(GNBG.o,1);
    case 2
        GNBG.ComponentSigma = GNBG.MinSigma + repmat(((GNBG.MaxSigma-GNBG.MinSigma)*rand),GNBG.o,1);
    case 3
        GNBG.ComponentSigma = [-1000;-950];%USER-DEFINED==>Number of rows must be equall to the number of components GNBG.o and number of columns is 1
    otherwise
        warning('Wrong number is chosen for GNBG.Overlapping.')
end
%% Defining the elements of diagonal elements of H for components (Configuring condition number)
GNBG.H_pattern       = 3;% (1) Condition number is 1 and all elements of principal diagonal of H are set to a user defined value GNBG.H_value
                         % (2) Condition number is 1 for all components but the elements of principal diagonal of H are different from a component to another and are randomly generated with uniform distribution within the range [GNBG.Lb_H, GNBG.Ub_H].
                         % (3) Condition number is random for all components the values of principal diagonal of the matrix H for each component are generated randomly within the range [GNBG.Lb_H, GNBG.Ub_H] using a uniform distribution.
                         % (4) Condition number is GNBG.Ub_H/GNBG.Lb_H for all components where two randomly selected elements on the principal diagonal of the matrix H are explicitly set to GNBG.Lb_H and GNBG.Ub_H. The remaining diagonal elements are generated randomly within the range [GNBG.Lb_H, GNBG.Ub_H]. These values follow a Beta distribution characterized by user-defined parameters alpha and beta, where 0 < alpha = beta <= 1.
                         % (5) Condition number is GNBG.Ub_H/GNBG.Lb_H for all components where a vector with GNBG.Dimension equally spaced values between GNBG.Lb_H and GNBG.Ub_H is generated. The linspace function is used to create a linearly spaced vector that includes both the minimum and maximum values. For each component, a randomly permutation of this vector is used. 
switch GNBG.H_pattern
    case 1
        GNBG.H_value = 1;
        GNBG.Component_H = GNBG.H_value .* ones(GNBG.o,GNBG.Dimension);
    case 2
        GNBG.Lb_H    = 1;
        GNBG.Ub_H    = 5;
        GNBG.Component_H   =  (GNBG.Lb_H + ((GNBG.Ub_H-GNBG.Lb_H)*rand(GNBG.o,1))).* ones(GNBG.o,GNBG.Dimension);
    case 3
        GNBG.Lb_H    = 1;
        GNBG.Ub_H    = 10^5;
        GNBG.Component_H   =  GNBG.Lb_H + ((GNBG.Ub_H-GNBG.Lb_H)*rand(GNBG.o,GNBG.Dimension));
    case 4
        GNBG.Lb_H    = 1;
        GNBG.Ub_H    = 10^5;
        GNBG.alpha   = 0.2;%Parameter for Beta distibution. With alpha = 1 and beta = 1, the beta distribution simplifies to a uniform distribution over the interval [0, 1].
        GNBG.beta    = GNBG.alpha;%alpha=beta because we need symmetric distribution
        GNBG.Component_H = GNBG.Lb_H + ((GNBG.Ub_H-GNBG.Lb_H)*betarnd(GNBG.alpha, GNBG.beta, [GNBG.o,GNBG.Dimension]));
        for ii = 1:GNBG.o
            random_indices = randperm(GNBG.Dimension, 2);
            GNBG.Component_H(ii, random_indices(1)) = GNBG.Lb_H;
            GNBG.Component_H(ii, random_indices(2)) = GNBG.Ub_H;
        end
    case 5
        GNBG.Lb_H    = 1;
        GNBG.Ub_H    = 10^5;
        GNBG.H_Values = linspace(GNBG.Lb_H, GNBG.Ub_H, GNBG.Dimension);
        for ii = 1:GNBG.o
            GNBG.Component_H(ii,:) = GNBG.H_Values(randperm(GNBG.Dimension));
        end
    otherwise
        warning('Wrong number is chosen for GNBG.H_pattern.')
end
%% Defining the parameters used in transformation T (Configuring modality of components and their basin local optima)
GNBG.MinMu            = 0.2;
GNBG.MaxMu            = 0.5;
GNBG.MinOmega         = 5;
GNBG.MaxOmega         = 50;
GNBG.localModalitySymmetry =3;%1 Unimodal, smooth, and regular components
                              %2 Multimodal symmetric components
                              %3 Multimodal asymmetric components
                              %4 Manually defined values
switch GNBG.localModalitySymmetry
    case 1
        GNBG.Mu    = zeros(GNBG.o,2);
        GNBG.Omega = zeros(GNBG.o,4);
    case 2
        GNBG.Mu    = repmat((GNBG.MinMu + (GNBG.MaxMu-GNBG.MinMu)*rand(GNBG.o,1)),1,2);
        GNBG.Omega = repmat((GNBG.MinOmega + (GNBG.MaxOmega-GNBG.MinOmega)*rand(GNBG.o,1)),1,4);
    case 3
        GNBG.Mu    = GNBG.MinMu + (GNBG.MaxMu-GNBG.MinMu)*rand(GNBG.o,2);
        GNBG.Omega = GNBG.MinOmega + (GNBG.MaxOmega-GNBG.MinOmega)*rand(GNBG.o,4);
    case 4
        GNBG.Mu    = [1,1];%USER-DEFINED==>Number of rows must be equall to the number of components GNBG.o and number of columns is 2
        GNBG.Omega = [10,10,10,10];%USER-DEFINED==>Number of rows must be equall to the number of components GNBG.o and number of columns is 4
    otherwise
        warning('Wrong number is chosen for GNBG.localModalitySymmetry.')
end
%% Defining the linearity of the basin of components
GNBG.MaxLambda           = 1;
GNBG.MinLambda           = 1;
GNBG.LambdaConfigMethod  = 1;%1 All lambda are set to GNBG.LambdaValue4ALL
                             %2 Randomly set lambda of each component in [GNBG.MinLambda,GNBG.MaxLambda]. Note that large ranges may result in existence of invisible components
switch GNBG.LambdaConfigMethod
    case 1
        GNBG.LambdaValue4ALL = 0.25;% GNBG.LambdaValue4ALL=0.5 for linear basin
                                 % GNBG.LambdaValue4ALL>0.5 for super-linear basin
                                 % 0<GNBG.LambdaValue4ALL<0.5 for sub-linear basin        
        GNBG.lambda  = GNBG.LambdaValue4ALL*ones(GNBG.o,1);
    case 2
        GNBG.lambda  = GNBG.MinLambda + (GNBG.MaxLambda-GNBG.MinLambda)*rand(GNBG.o,1);
    otherwise
        warning('Wrong number is chosen for GNBG.LambdaConfigMethod.')
end
%% Defining variable interaction structure and generating rotation matrices
GNBG.MinAngle = -pi;
GNBG.MaxAngle = pi;
GNBG.Rotation = 3;%(1) Without rotation
                  %(2) Random Rotation for all components==> Fully Connected Interaction with random angles for each plane of each component
                  %(3) Random Rotation for all components==>For each component, interactions are defined randomly based on connection probability threshold GNBG.ConnectionProbability and random angles. GNBG.ConnectionProbability can be different for each component.
                  %(4) Rotation based on the random Angle for each component (fully connected with an angle for all planes in each component)
                  %(5) Rotation based on the specified Angle for all components (fully connected with an angle for all planes in all component)
                  %(6) Rotation based on the random Angles to generate chain-like variable interaction structure
                  %(7) Generating partialy separable variable interactiion structure with user defined sizes and angles for each group of variables
switch GNBG.Rotation
    case 1
        GNBG.RotationMatrix = NaN(GNBG.Dimension,GNBG.Dimension,GNBG.o);
        for ii=1 : GNBG.o
            GNBG.RotationMatrix(:,:,ii) = eye(GNBG.Dimension);
        end
    case 2
        GNBG.RotationMatrix = NaN(GNBG.Dimension,GNBG.Dimension,GNBG.o);
        for ii=1 : GNBG.o
            GNBG.ThetaMatrix = zeros(GNBG.Dimension);
            upperTriangle = triu(true(GNBG.Dimension), 1); 
            GNBG.ThetaMatrix(upperTriangle) = GNBG.MinAngle + (GNBG.MaxAngle - GNBG.MinAngle) * rand(sum(upperTriangle(:)), 1);
            [GNBG.RotationMatrix(:,:,ii)] = Rotation(GNBG.ThetaMatrix,GNBG.Dimension);
        end
    case 3
        GNBG.RotationMatrix = NaN(GNBG.Dimension,GNBG.Dimension,GNBG.o);
        GNBG.MinConProb     = 0.75;%If an equal value of GNBG.ConnectionProbability for all components is needed, then set GNBG.MinConProb and GNBG.MaxConProb to an identical value.
        GNBG.MaxConProb     = 0.75;
        GNBG.ConnectionProbability = GNBG.MinConProb + (GNBG.MaxConProb-GNBG.MinConProb)*rand(GNBG.o,1);
        for ii=1 : GNBG.o 
            GNBG.ThetaMatrix=zeros(GNBG.Dimension);
            upperTriangle = triu(rand(GNBG.Dimension)<GNBG.ConnectionProbability(ii),1); 
            GNBG.ThetaMatrix(upperTriangle) = GNBG.MinAngle + (GNBG.MaxAngle - GNBG.MinAngle) * rand(sum(upperTriangle(:)), 1);
            [GNBG.RotationMatrix(:,:,ii)] = Rotation(GNBG.ThetaMatrix,GNBG.Dimension);
        end    
    case 4
        GNBG.RotationMatrix = NaN(GNBG.Dimension,GNBG.Dimension,GNBG.o);
        GNBG.Angle = GNBG.MinAngle + (GNBG.MaxAngle-GNBG.MinAngle)*rand(GNBG.o,1);
        for ii=1 : GNBG.o
            GNBG.ThetaMatrix=zeros(GNBG.Dimension);
            upperTriangle = triu(true(GNBG.Dimension), 1); % Logical matrix for upper triangle
            GNBG.ThetaMatrix(upperTriangle) = GNBG.Angle(ii);
            [GNBG.RotationMatrix(:,:,ii)] = Rotation(GNBG.ThetaMatrix,GNBG.Dimension);
        end
    case 5
        GNBG.RotationMatrix = NaN(GNBG.Dimension,GNBG.Dimension,GNBG.o);
        GNBG.Angle = (1.48353)*ones(GNBG.o,1);%User specified angle for all planes in all components
        for ii=1 : GNBG.o
            GNBG.ThetaMatrix=zeros(GNBG.Dimension);
            upperTriangle = triu(true(GNBG.Dimension), 1); % Logical matrix for upper triangle
            GNBG.ThetaMatrix(upperTriangle) = GNBG.Angle(ii);
            [GNBG.RotationMatrix(:,:,ii)] = Rotation(GNBG.ThetaMatrix,GNBG.Dimension);
        end
    case 6
        GNBG.RotationMatrix = NaN(GNBG.Dimension,GNBG.Dimension,GNBG.o);
        for ii=1 : GNBG.o
            GNBG.ThetaMatrix = zeros(GNBG.Dimension);
            GNBG.ThetaMatrix(sub2ind(size(GNBG.ThetaMatrix), 1:GNBG.Dimension-1, 2:GNBG.Dimension)) = GNBG.MinAngle + (GNBG.MaxAngle - GNBG.MinAngle) .* rand(GNBG.Dimension-1, 1);
            [GNBG.RotationMatrix(:,:,ii)] = Rotation(GNBG.ThetaMatrix,GNBG.Dimension);
        end
    case 7
        GNBG.RotationMatrix = NaN(GNBG.Dimension,GNBG.Dimension,GNBG.o);
        S=[10,10,10];%USER-DEFINED==>The size of S shows the number of groups of variables that have interactions among themselves and the values of S show the size of each group. Sum of elements in S must be less than or equal to GNBG.Dimension. 
        Theta=[pi/4,3*pi/4,pi/8];%USER-DEFINED==>Shows the angle of rotation used for each group of variables in S. Size of Theta must be the same as S.
        if sum(S) > GNBG.Dimension || length(S) ~= length(Theta)
            error('The sum of elements in S exceeds GNBG.Dimension, or the size of S and Theta are not equal.');
        end
        for ii=1 : GNBG.o
            GNBG.ThetaMatrix = zeros(GNBG.Dimension);
            allVars = randperm(GNBG.Dimension);
            groupStart = 1;
            for jj = 1:length(S)
                groupEnd = groupStart + S(jj) - 1;
                groupVars = allVars(groupStart:groupEnd);
                for var1 = groupVars
                    for var2 = groupVars
                        if var1 < var2 % Only for elements above the diagonal
                            GNBG.ThetaMatrix(var1, var2) = Theta(jj);%Can be replaced by a random angle by replacing Theta(i) with GNBG.MinAngle+(GNBG.MaxAngle-GNBG.MinAngle)*rand
                        end
                    end
                end
                groupStart = groupEnd + 1;
            end
            [GNBG.RotationMatrix(:,:,ii)] = Rotation(GNBG.ThetaMatrix,GNBG.Dimension);
        end
    otherwise
        warning('Wrong number is chosen for GNBG.Rotation.')
end
%% Parameters used for storing the results
[GNBG.OptimumValue,GlobalOptimumID] = min(GNBG.ComponentSigma);
GNBG.OptimumPosition  = GNBG.Component_MinimumPosition(GlobalOptimumID,:);
GNBG.FEhistory        = inf(1,GNBG.MaxEvals);
GNBG.FE               = 0;
GNBG.AcceptanceReachPoint = inf;
GNBG.BestFoundResult      = inf;
GNBG.BestAtFirstLine      = NaN;
GNBG.BestAtSecondLine     = NaN;
end
%% Rotation matrix generator function
function R = Rotation(teta,Dimension)
R = eye(Dimension);
for p=1 : (Dimension-1)
    for q=(p+1) : (Dimension)
        if teta(p,q)~=0
            G = eye(Dimension);
            G(p,p) = cos(teta(p,q));
            G(q,q) = cos(teta(p,q));
            G(p,q) = -sin(teta(p,q));
            G(q,p) = sin(teta(p,q));
            R = R*G;
        end
    end
end
end