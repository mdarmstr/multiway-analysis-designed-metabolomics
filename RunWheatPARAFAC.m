%% PARAFAC Example on Wheat Data
%
% Parallel Factor Analysis (PARAFAC) example with the data collected in Warth, B. et al. (2014).
% Metabolomics, 11(3), 722�738. Data were downloaded from the MetaboLights
% metabolomics public data repository (www.ebi.ac.uk/\-metabolights, with
% accession number MTBLS112).  Experiments aimed at identifying changes in
% the metabolome of wheat (Triticum aestivum) induced by deoxynivalenol
% (DON), a mycotxin produced by the infestant Fusarium graminearum and
% related species causing the devastating plant disease Fusarium head
% blight. In the study, four wheat genotypes with known varying resistance
% to Fusarium were treated with either DON or water control and harvested
% after 0, 12, 24, 48 and 96 hours after treatment. Target GC-MS profiling
% was used to quantify an array of 57 metabolites. The resulting data
% matrix X has dimensions 296 x 57.
%
% coded by: Michael Sorochan Armstrong (mdarmstr@ugr.es), and Jose Camacho
% Paez (jcamachop@ugr.es)
% last modification: 22/Dec/2022
%
% Copyright (C) 2023  University of Granada, Granada
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%Please run RunWheatUPCA first to acquire the necessary data for this script.

addpath 'MEDA\MEDA-Toolbox-master\'
addpath 'nway320 exchange'\

clear all %#ok
close all
load wheat

X_m = X - mean(X);


%% As per Pepe's suggestion, and working towards a general formulation, let's create a tensor using a for loop
utime = unique(time); % factor time
for i=1:length(time)
    for j=1:length(utime)
        if strcmp(time{i},utime(j))
            ytim(i)=j;
        end
    end
end

utrait = unique(trait); % factor trait
for i=1:length(trait)
    for j=1:length(utrait)
        if strcmp(trait{i},utrait(j))
            ytra(i)=j;
        end
    end
end

utreat = unique(treatment); % factor treatment
for i=1:length(treatment)
    for j=1:length(utreat)
        if strcmp(treatment{i},utreat(j))
            ytre(i)=j;
        end
    end
end

urep = unique(replicate); % factor replicate
for ii = 1:length(replicate)
    for jj = 1:length(urep)
        if strcmp(replicate{ii},urep{jj})
            yrep(ii) = jj;
        end
    end
end

F = [ytre',yrep',ytim',ytra'];

%% Create the tensor in a for loop.

tnsr = zeros(2,5,5,4,58); %Treatment, Replicate, Time, Trait, Metabolite

for ii = 1:length(utreat)
    for jj = 1:length(urep)
        for kk = 1:length(utime)
            for ll = 1:length(utrait)
                for mm = 1:length(var_l)
                    tnsr(ii,jj,kk,ll,mm) = X_m(F(:,1)==ii & F(:,2) == jj & F(:,3) == kk & F(:,4) == ll,mm);
                end
            end
        end
    end
end

X = tnsr;

%% What is the percent variance explained by the PARAFAC model, including the replicates?
pftest(3,X,5,[0 0 0 0 NaN]);
[~,~,err,corr] = parafac(X,1,[0,0,0,2,0,0]);
disp((sum(X(:).^2) - err)/ sum(X(:).^2) * 100) %32.97
disp(corr) %100

pause()

X_cell = cell(1,size(urep,1));
Fstrct = cell(1,size(urep,1));

for ii = 1:size(urep,1)
    X_cell{ii} = squeeze(X(:,ii,:,:,:));
end

%Test for number of components - 1 looks good.
[ssX,Corco] = pftest(3,X,5,[0 0 0 0 NaN]);
saveas(gcf,'Fig/parafactest.eps','epsc');
pvar_parafac = zeros(5,1);
corr_parafac = zeros(5,1);

for ii = 1:size(urep,1)
    %Build a 1-component PARAFAC model without any additional scaling and no constraints.
    [F,iter,err,corr] = parafac(X_cell{1,ii},1,[0,0,0,2,0,0]);
    pvar_parafac(ii) = (sum(tnsr(:).^2) - err)/ sum(tnsr(:).^2) * 100;
    corr_parafac(ii) = corr;
    Fstrct{1,ii} = F;
end

%Plotting the treatment mode

for ii = 1:size(urep,1)
    lds(ii,:) = Fstrct{ii}{1};
end
md_mean = mean(lds);
md_stdv = std(lds);

plot_vec(md_mean,utreat,utreat); hold on;
errorbar(md_mean,md_stdv,'.','CapSize',18,'Color','k'); hold off;
title('PARAFAC loadings of treatment mode')
saveas(gcf,'Fig/scores_parafac.eps','epsc');

clear lds

%plotting the time mode
for ii = 1:size(urep,1)
    lds(ii,:) = Fstrct{ii}{2};
end
md_mean = mean(lds);
md_stdv = std(lds);

plot_vec(md_mean,utime); hold on;
errorbar(md_mean,md_stdv,'.','CapSize',18,'Color','k'); hold off;
title('PARAFAC loadings of time mode')
saveas(gcf,'Fig/loadings_parafactime.eps','epsc');

clear lds

%plotting the trait mode
for ii = 1:size(urep,1)
    lds(ii,:) = Fstrct{ii}{3};
end
md_mean = mean(lds);
md_stdv = std(lds);

%Plotting the trait mode
plot_vec(md_mean,utrait); hold on;
errorbar(md_mean,md_stdv,'.','CapSize',18,'Color','k'); hold off;
title('PARAFAC loadings of trait mode')
saveas(gcf,'Fig/loadings_parafactrait.eps','epsc');

clear lds

%plotting the trait mode
for ii = 1:size(urep,1)
    lds(ii,:) = Fstrct{ii}{4};
end
md_mean = mean(lds);
md_stdv = std(lds);

%Plotting the metabolite mode
plot_vec(md_mean,var_l); hold on;
errorbar(md_mean,md_stdv,'.','CapSize',6,'Color','k'); hold off;
title('PARAFAC loadings of metabolite mode')
saveas(gcf,'Fig/loadings_parafacmet.eps','epsc');

%% Plotting the error as a function of the different replicates
plot_vec([corr_parafac, pvar_parafac],{1,2,3,4,5});
title('PARAFAC Model Performance per Replicate')
xlabel('Replicates');
ylabel('% VAR EXP / CONCORDIA')
legend({'CONCORDIA', 'VAR. EXPLAINED'},'Location','southeast')
saveas(gcf,'Fig/parafac_performance.eps','epsc');

clear lds

save("parafac_met.mat","md_mean");

