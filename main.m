%***********************************GNBG+DE***********************************
%Author: Danial Yazdani
%Last Edited: December 14, 2023
%Title: Generalized Numerical Benchmark Generator
% This MATLAB code implements the Differential Evolution (DE) optimization 
% algorithm (DE/rand/1/bin) to solve problem instances generated
% by the Generalized Numerical Benchmark Generator (GNBG). The GNBG 
% parameter setting for each problem instance is configured in 'BenchmarkGenerator.m'. 
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
close all;clear all;clc; %#ok<CLALL> 
RunNumber = 3;
Error = NaN(3,RunNumber);
AcceptancePoints = NaN(1,RunNumber);
for RunCounter=1 : RunNumber
    disp(['RunCounter=', num2str(RunCounter)]);
    %% Preparation and loading of the GNBG parameters based on the chosen problem instance
    clear GNBG DE;
    rng(1);% GNBG uses multiple random numbers to generate problem instances which affect some problem characteristics. An identical random seed should be used for all runs to generate identical problem instances for each experiment.
    GNBG = BenchmarkGenerator;
    rng('shuffle');%Set a random seed for the optimizer
    %% Optimizer part
    DE.PopulationSize = 100;
    DE.Dimension = GNBG.Dimension;
    DE.LB = GNBG.MinCoordinate;
    DE.UB = GNBG.MaxCoordinate;
    DE.X = DE.LB + (DE.UB - DE.LB)*rand(DE.PopulationSize,DE.Dimension);
    [DE.FitnessValue,GNBG] = fitness(DE.X,GNBG);
    DE.Donor = NaN(DE.PopulationSize,DE.Dimension);
    DE.Cr = 0.9;
    DE.F = 0.5;
    [~,DE.BestID] = min(DE.FitnessValue);
    DE.BestPosition = DE.X(DE.BestID,:);
    DE.BestValue = DE.FitnessValue(DE.BestID);
    %% main loop of the optimizer
    while 1
        [~,DE.BestID] = min(DE.FitnessValue);
        if DE.FitnessValue(DE.BestID)<DE.BestValue
            DE.BestPosition = DE.X(DE.BestID,:);
            DE.BestValue = DE.FitnessValue(DE.BestID);
        end
        %% Mutation
        R = NaN(DE.PopulationSize,3);
        for ii=1 : DE.PopulationSize
            tmp = randperm(DE.PopulationSize);
            tmp(tmp==ii)=[];
            R(ii,:) = tmp(1:3);
        end
        DE.Donor = DE.X(R(:,1),:) + DE.F.*(DE.X(R(:,2),:)-DE.X(R(:,3),:));%DE/rand/1
        %% Crossover==>binomial
        DE.OffspringPosition = DE.X;%U
        K = sub2ind([DE.PopulationSize,DE.Dimension],(1:DE.PopulationSize)',randi(DE.Dimension,[DE.PopulationSize,1]));
        DE.OffspringPosition(K) = DE.Donor(K);
        CrossoverBinomial = rand(DE.PopulationSize,DE.Dimension)<repmat(DE.Cr,DE.PopulationSize,DE.Dimension);
        DE.OffspringPosition(CrossoverBinomial) = DE.Donor(CrossoverBinomial);
        %% boundary checking
        LB_tmp1 = DE.OffspringPosition<DE.LB;
        LB_tmp2 = ((DE.LB + DE.X).*LB_tmp1)/2;
        DE.OffspringPosition(LB_tmp1) = LB_tmp2(LB_tmp1);
        UB_tmp1 = DE.OffspringPosition>DE.UB;
        UB_tmp2 = ((DE.UB + DE.X).*UB_tmp1)/2;
        DE.OffspringPosition(UB_tmp1) = UB_tmp2(UB_tmp1);
        [DE.OffspringFitness, GNBG] = fitness(DE.OffspringPosition(:,1:DE.Dimension), GNBG);
        %% Selection==>greedy
        better = DE.OffspringFitness < DE.FitnessValue;
        DE.X(better,:) = DE.OffspringPosition(better,:);
        DE.FitnessValue(better) = DE.OffspringFitness(better);
        
        if  GNBG.FE >= GNBG.MaxEvals%When termination criteria has been met
            break;
        end
    end
    %% Storing results of each run
    Error(1,RunCounter) = abs(GNBG.BestAtFirstLine - GNBG.OptimumValue);
    Error(2,RunCounter) = abs(GNBG.BestAtSecondLine - GNBG.OptimumValue);
    Error(3,RunCounter) = abs(GNBG.BestFoundResult - GNBG.OptimumValue);
    for jj=2:GNBG.MaxEvals
       if GNBG.FEhistory(jj-1)<GNBG.FEhistory(jj)
           GNBG.FEhistory(jj) = GNBG.FEhistory(jj-1);
       end
    end
    ConvBhv(RunCounter,:) = abs(GNBG.FEhistory(1:GNBG.MaxEvals) - GNBG.OptimumValue); %#ok<SAGROW>
    AcceptancePoints(RunCounter) = GNBG.AcceptanceReachPoint;
end
%% Output
nonInfIndices = isfinite(AcceptancePoints);
nonInfValues = AcceptancePoints(nonInfIndices);
disp(['Average FE to reach acceptance result: ', num2str(mean(nonInfValues)),'(',num2str(std(nonInfValues)),')']);
disp(['Acceptance Ratio: ', num2str((sum(nonInfIndices) / length(AcceptancePoints)) * 100)]);
disp(['FirstPointError: ', num2str(mean(Error(1,:))),'(',num2str(std(Error(1,:))),')']);
disp(['SecondPointError: ', num2str(mean(Error(2,:))),'(',num2str(std(Error(2,:))),')']);
disp(['LastPointError: ', num2str(mean(Error(3,:))),'(',num2str(std(Error(3,:))),')']);

%% Convergence plot generation
% MeanConvBhv = mean(ConvBhv);
% interval = GNBG.MaxEvals / 10; 
% numbers = linspace(interval, GNBG.MaxEvals - interval, 10);
% f=semilogy(MeanConvBhv,'k', 'LineWidth',1, 'Marker', 'd', 'MarkerIndices',round(numbers),'DisplayName','Main optima # = ');
% xlabel('Fitness Evaluation Number')
% ylabel('Error')
% ylim([0,MeanConvBhv(1)])
% grid on
% hold on 
% line([0, GNBG.MaxEvals], [GNBG.AcceptanceThreshold, GNBG.AcceptanceThreshold], 'Color', 'red', 'LineStyle', '--', 'HandleVisibility', 'off', 'LineWidth',1)
% xline_value = GNBG.FirstPoint;
% xline(xline_value, 'Color', 'blue', 'LineStyle', '--', 'HandleVisibility', 'off', 'LineWidth',1)
% xline_value = GNBG.SecondPoint;
% xline(xline_value, 'Color', 'blue', 'LineStyle', '--', 'HandleVisibility', 'off', 'LineWidth',1)
% set(gcf,'OuterPosition',[150 150 600 550]);
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% xline_value = round(0.75*GNBG.MaxEvals);
% xline(xline_value, 'Color', 'blue', 'LineStyle', '--')
% formatOut = 'mm-dd-yy-HH-MM-SS';
% Figruename = strcat('DE-',datestr(now,formatOut));
% saveas(gcf,Figruename,'epsc')
% saveas(gcf,strcat(Figruename,'.jpg'))
% saveas(gcf,strcat(Figruename,'.fig'))
% disp([num2str(mean(Error(1,:))),'(',num2str(std(Error(1,:))),') & ', num2str(mean(Error(2,:))),'(',num2str(std(Error(2,:))),') & ', num2str(mean(Error(3,:))),'(',num2str(std(Error(3,:))),') & ', num2str(mean(nonInfValues)),'(',num2str(std(nonInfValues)),') & ', num2str((sum(nonInfIndices) / length(AcceptancePoints)) * 100), ' \\']);
