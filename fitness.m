%***********************************GNBG***********************************
%Author: Danial Yazdani
%Last Edited: December 14, 2023
%Title: Generalized Numerical Benchmark Generator
% --------
%Description: 
%          This function includes the objective function of GNBG based on 
%          parameter settings defined by the user and stored in 'GNBG' structure.
%          In addition, the results of the algorithms is gathered in this function
%          and stored in 'GNBG' structure. 
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
function [result,GNBG] = fitness(X,GNBG)
[SolutionNumber,~] = size(X);
result = NaN(SolutionNumber,1);
for jj=1 : SolutionNumber
    x = X(jj,:)';
    f=NaN(1,GNBG.o);
    for k=1 : GNBG.o
        a = Transform((x - GNBG.Component_MinimumPosition(k,:)')'*GNBG.RotationMatrix(:,:,k)',GNBG.Mu(k,:),GNBG.Omega(k,:));
        b = Transform(GNBG.RotationMatrix(:,:,k) * (x - GNBG.Component_MinimumPosition(k,:)'),GNBG.Mu(k,:),GNBG.Omega(k,:));
        f(k) = GNBG.ComponentSigma(k) + ( a * diag(GNBG.Component_H(k,:)) * b)^GNBG.lambda(k);
    end
    result(jj) = min(f);
    if GNBG.FE > GNBG.MaxEvals
        return;
    end
    GNBG.FE = GNBG.FE + 1;
    GNBG.FEhistory(GNBG.FE) = result(jj);
    %%
    if GNBG.BestFoundResult > result(jj)
        GNBG.BestFoundResult = result(jj);
    end
    if (abs(GNBG.FEhistory(GNBG.FE) - GNBG.OptimumValue)) < GNBG.AcceptanceThreshold && isinf(GNBG.AcceptanceReachPoint)
        GNBG.AcceptanceReachPoint = GNBG.FE;
    end
    if GNBG.FE == GNBG.FirstPoint
        GNBG.BestAtFirstLine = min(GNBG.FEhistory(1:GNBG.FE));
    end
    if GNBG.FE == GNBG.SecondPoint
        GNBG.BestAtSecondLine = min(GNBG.FEhistory(1:GNBG.FE));
    end
    %%
end
end

function Y = Transform(X,Alpha,Beta)
Y = X;
tmp = (X > 0);
Y(tmp) = log(X(tmp));
Y(tmp) = exp(Y(tmp) + Alpha(1)*(sin(Beta(1).*Y(tmp)) + sin(Beta(2).*Y(tmp))));
tmp = (X < 0);
Y(tmp) = log(-X(tmp));
Y(tmp) = -exp(Y(tmp) + Alpha(2)*(sin(Beta(3).*Y(tmp)) + sin(Beta(4).*Y(tmp))));
end