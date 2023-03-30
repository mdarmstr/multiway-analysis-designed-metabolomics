%% PARAFASCA Example on Wheat Data
%
% Parallel Factor Analysis (PARAFAC) example with the data collected in Warth, B. et al. (2014).
% Metabolomics, 11(3), 722-738. Data were downloaded from the MetaboLights
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



%% Add the path to the necessary toolboxes, if they have not already been added. Load the data.
addpath 'MEDA\MEDA-Toolbox-master'
addpath 'nway320 exchange'\

clear all %#ok
close all
load wheat

Xm = X - mean(X);

%% Building the design matrix

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

%% Assuming the replicates are not a significant experimental factor, build the design matrix without them
Fp = [ytre',ytim',ytra'];

%% Perform a generalised linear model of the data with the design matrix.
[tbl,paranovao] = parglm(Xm, Fp, [1,2; 1,3; 2,3],0); %Already mean-centred

disp(tbl)

pause(1)

%% Build a matrix from the significant factors and interactions.
Xr = paranovao.factors{1}.matrix + paranovao.factors{2}.matrix + paranovao.factors{3}.matrix ...
        + paranovao.interactions{1}.matrix + paranovao.interactions{2}.matrix + paranovao.interactions{3}.matrix;

Tr = zeros(2,5,5,4,58);
Tm = zeros(2,5,5,4,58);

for ii = 1:length(utreat)
    for jj = 1:length(urep)
        for kk = 1:length(utime)
            for ll = 1:length(utrait)
                for mm = 1:length(var_l)
                    Tr(ii,jj,kk,ll,mm) = Xr(F(:,1)==ii & F(:,2) == jj & F(:,3) == kk & F(:,4) == ll,mm);
                    Tm(ii,jj,kk,ll,mm) = Xm(F(:,1)==ii & F(:,2) == jj & F(:,3) == kk & F(:,4) == ll,mm);
                end
            end
        end
    end
end

%% Test for the number of PARAFAC components using one replicate

To = squeeze(Tr(:,ii,:,:,:));

[ssX,Corco] = pftest(3,To,5,[0,0,0,2,NaN]);
saveas(gcf,'Fig/parafascatest.eps','epsc');

%% Build a one-component PARAFAC model on the response tensor

[Fc,iter,err,corr] = parafac(To,1,[0,0,0,2,0,0]);

%Preallocating space for replicates.
Te = cell(5,1);
Fp_1 = cell(5,1);
Fp_2 = cell(5,1);
Fp_3 = cell(5,1);
Fp_4 = cell(5,1);


for ii = 1:5

    Te{ii} = To - squeeze(Tm(:,ii,:,:,:));

    %Getting the first treatment factor
    Z = krb(Fc{4},krb(Fc{3},Fc{2}));
    Xi_jkl = (nshape(To,1) + nshape(Te{ii},1));
    Fp_1{ii} = Xi_jkl * Z * pinv(Z'*Z);

    %Getting the time factor
    Z = krb(Fc{4},krb(Fc{3},Fc{1}));
    Xj_ikl = (nshape(To,2) + nshape(Te{ii},2));
    Fp_2{ii} = Xj_ikl * Z * pinv(Z'*Z);

    %Getting the trait factor
    Z = krb(Fc{4},krb(Fc{2},Fc{1}));
    Xk_ijl = (nshape(To,3) + nshape(Te{ii},3));
    Fp_3{ii} = Xk_ijl * Z * pinv(Z'*Z);

    %Getting the metabolite factor just in case it's relevant
    Z = krb(Fc{3},krb(Fc{2},Fc{1}));
    Xl_ijk = (nshape(To,4) + nshape(Te{ii},4));
    Fp_4{ii} = Xl_ijk * Z * pinv(Z'*Z);

end

%% Plotting results

%Plotting the treatment mode

for ii = 1:size(urep,1)
    lds1(ii,:) = Fp_1{ii}; %#ok it is a small loop
end
md_mean = mean(lds1);
md_stdv = std(lds1);

plot_vec(md_mean,utreat,utreat); hold on;
errorbar(md_mean,md_stdv,'.','CapSize',18,'Color','k'); hold off;
title('PARAFASCA loadings of treatment mode')
saveas(gcf,'Fig/scores_parafasca.eps','epsc');

%plotting the time mode
for ii = 1:size(urep,1)
    lds2(ii,:) = Fp_2{ii}; %#ok
end
md_mean = mean(lds2);
md_stdv = std(lds2);

plot_vec(md_mean,utime); hold on;
errorbar(md_mean,md_stdv,'.','CapSize',18,'Color','k'); hold off;
title('PARAFASCA loadings of time mode')
saveas(gcf,'Fig/loadings_parafascatime.eps','epsc');

%plotting the trait mode
for ii = 1:size(urep,1)
    lds3(ii,:) = Fp_3{ii}; %#ok
end
md_mean = mean(lds3);
md_stdv = std(lds3);

%Plotting the trait mode
plot_vec(md_mean,utrait); hold on;
errorbar(md_mean,md_stdv,'.','CapSize',18,'Color','k'); hold off;
title('PARAFASCA loadings of trait mode')
saveas(gcf,'Fig/loadings_parafascatrait.eps','epsc');

%plotting the trait mode
for ii = 1:size(urep,1)
    lds4(ii,:) = Fp_4{ii}; %#ok
end
md_mean = mean(lds4);
md_stdv = std(lds4);

%plotting the metabolite mode
plot_vec(md_mean,var_l); hold on;
errorbar(md_mean,md_stdv,'.','CapSize',6,'Color','k'); hold off;
title('PARAFASCA loadings of metabolite mode')
saveas(gcf,'Fig/loadings_parafascamet.eps','epsc');

%% Calculating on percent variance explained on one replicate
pvar = (sum(To(:).^2) - err)/ sum(To(:).^2) * 100;
disp(pvar)
disp(corr)

save("parafasca_met.mat","md_mean")
 
