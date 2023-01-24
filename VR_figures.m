%% Copepod Variance Ratio Calculations
% Isabel Honda
% January 2023

%% Community Synchrony - GAM Generated Timeseries

clc; clear; close all

locs = ["MVCO","L1","L2","L3","L4","L5","L6","L7","L8","L9","L10","L11"];
comm_VRs_entire = zeros(length(locs),2);
comm_VRs_int = zeros(length(locs),2);
comm_VRs_seas = zeros(length(locs),2);

figure(1)
clf
set(gcf,'color','w');

for i=1:length(locs)
    calfin = readtable(sprintf('Transect_timeseries/calfin_%s.csv',locs(i)));
    pseudo = readtable(sprintf('Transect_timeseries/pseudo_%s.csv',locs(i)));
    ctyp = readtable(sprintf('Transect_timeseries/ctyp_%s.csv',locs(i)));
    comm_loc_entire = [calfin.log10_calfinGAMest pseudo.log10_pseudoGAMest ctyp.log10_ctypGAMest];
    comm_VRs_entire(i) = variance_ratio(comm_loc_entire);
    [lb,ub] = VR_sigTest(comm_loc_entire);
    if comm_VRs_entire(i,1) > ub || comm_VRs_entire(i,1) < lb
        comm_VRs_entire(i,2) = 1;
    end

    avgYear_calfin = accumarray(calfin.year,calfin.log10_calfinGAMest,[],@mean);
    avgYear_pseudo = accumarray(pseudo.year,pseudo.log10_pseudoGAMest,[],@mean);
    avgYear_ctyp = accumarray(ctyp.year,ctyp.log10_ctypGAMest,[],@mean);
    comm_loc_int = [nonzeros(avgYear_calfin) nonzeros(avgYear_pseudo) nonzeros(avgYear_ctyp)];
    comm_VRs_int(i) = variance_ratio(comm_loc_int);
    [lb,ub] = VR_sigTest(comm_loc_int);
    if comm_VRs_int(i,1) > ub || comm_VRs_int(i,1) < lb
        comm_VRs_int(i,2) = 1;
    end

    seas_calfin = accumarray(calfin.julian,calfin.log10_calfinGAMest,[],@mean);
    seas_pseudo = accumarray(pseudo.julian,pseudo.log10_pseudoGAMest,[],@mean);
    seas_ctyp = accumarray(ctyp.julian,ctyp.log10_ctypGAMest,[],@mean);
    comm_loc_seas = [seas_calfin seas_pseudo seas_ctyp];
    comm_VRs_seas(i) = variance_ratio(comm_loc_seas);
    [lb,ub] = VR_sigTest(comm_loc_seas);
    if comm_VRs_seas(i,1) > ub || comm_VRs_seas(i,1) < lb
        comm_VRs_seas(i,2) = 1;
    end

    if i==1
        subplot(3,3,1)
        plot(calfin.year + calfin.julian/365,comm_loc_entire,'LineWidth',2)
        xlim([1977 2019])
        legend('Calfin','Pseudo','Ctyp')
        title('Entire Timeseries at MVCO')
        ylabel('log_{10}(N+1) m^{-3}')
        xlabel('Year')

        subplot(3,3,4)
        plot(1977:2019,comm_loc_int,'LineWidth',2)
        xlim([1977 2019])
        legend('Calfin','Pseudo','Ctyp')
        title('Interannual Variability at MVCO')
        ylabel('log_{10}(N+1) m^{-3}')
        xlabel('Year')
        
        subplot(3,3,7)
        plot(1:365,comm_loc_seas,'LineWidth',2)
        xlim([1 365])
        legend('Calfin','Pseudo','Ctyp')
        title('Seasonal Variability at MVCO')
        ylabel('log_{10}(N+1) m^{-3}')
        xlabel('Day of Year')

    elseif i==12
        subplot(3,3,2)
        plot(calfin.year + calfin.julian/365,comm_loc_entire,'LineWidth',2)
        xlim([1977 2019])
        legend('Calfin','Pseudo','Ctyp')
        title('Entire Timeseries at L11')
        ylabel('log_{10}(N+1) m^{-3}')
        xlabel('Year')

        subplot(3,3,5)
        plot(1977:2019,comm_loc_int,'LineWidth',2)
        xlim([1977 2019])
        legend('Calfin','Pseudo','Ctyp')
        title('Interannual Variability at L11')
        ylabel('log_{10}(N+1) m^{-3}')
        xlabel('Year')
        
        subplot(3,3,8)
        plot(1:365,comm_loc_seas,'LineWidth',2)
        xlim([1 365])
        legend('Calfin','Pseudo','Ctyp')
        title('Seasonal Variability at L11')
        ylabel('log_{10}(N+1) m^{-3}')
        xlabel('Day of Year')
    end
end


subplot(3,3,3)
plot(1:length(locs),comm_VRs_entire(:,1),'LineWidth',2,'Color','#a87d60')
hold on
for i=1:length(locs)
    if comm_VRs_entire(i,2) == 1
        h(1)=plot(i,comm_VRs_entire(i,1),'*','Color','#a87d60','Linewidth',2)
    end   
end
set(gca, 'Xtick',1:length(locs),'XTickLabel',locs);
xlabel('Station')
ylabel('Variance Ratio')
title('Entire Timeseries VR')
xlim([1 length(locs)])
legend(h,'Statistically Significant (95% CI)','Location','southeast')

subplot(3,3,6)
plot(1:length(locs),comm_VRs_int(:,1),'LineWidth',2,'Color','#a87d60')
hold on
for i=1:length(locs)
    if comm_VRs_int(i,2) == 1
        h(1)=plot(i,comm_VRs_int(i,1),'*','Color','#a87d60','Linewidth',2)
    end   
end
set(gca, 'Xtick',1:length(locs),'XTickLabel',locs);
xlabel('Station')
ylabel('Variance Ratio')
title('Interannual Variability VR')
xlim([1 length(locs)])
legend(h,'Statistically Significant (95% CI)','Location','southeast')

subplot(3,3,9)
plot(1:length(locs),comm_VRs_seas(:,1),'LineWidth',2,'Color','#a87d60')
hold on
for i=1:length(locs)
    if comm_VRs_seas(i,2) == 1
        h(1)=plot(i,comm_VRs_seas(i,1),'*','Color','#a87d60','Linewidth',2)
    end   
end
set(gca, 'Xtick',1:length(locs),'XTickLabel',locs);
xlabel('Station')
ylabel('Variance Ratio')
title('Seasonal Variability VR')
xlim([1 length(locs)])
legend(h,'Statistically Significant (95% CI)','Location','southeast')


%% Community Synchrony - Raw EcoMon data (no GAM)


clc; clear; close all
s22 = readtable('EcoMon_Agg/strata_22.csv');
s23 = readtable('EcoMon_Agg/strata_23.csv');
s24 = readtable('EcoMon_Agg/strata_24.csv');
s25 = readtable('EcoMon_Agg/strata_25.csv');
oldStrats = load('NES_polygons46.mat');

strat = {s25, s24, s23, s22};

VRs = zeros(length(strat),2);

for i=1:length(strat)
    VRs(i,1) = variance_ratio([strat{i}.calfin_100m3 strat{i}.pseudo_100m3, strat{i}.ctyp_100m3]);
    [lb,ub] = VR_sigTest([strat{i}.calfin_100m3 strat{i}.pseudo_100m3, strat{i}.ctyp_100m3]);
    if VRs(i,1) > ub || VRs(i,1) < lb
        VRs(i,2) = 1;
    end
end


figure(2)
clf
set(gcf,'color','w');
subplot(4,3,[3 6])
for i=22:25
    plot(oldStrats.shape(i),'EdgeColor','black','FaceColor','white');  
    hold on
    [xn,yn] = centroid(oldStrats.shape(i));
    if i==24
        xn = xn+0.2;
    end
    text(xn,yn-0.1,num2str(i));
end
xlabel('Longitude')
ylabel('Latitude')
title('Strata')

titles = {'Strata 25 - Coastal','Strata 24', 'Strata 23', 'Strata 22 - Offshore'};

i=1;
for j=1:3:10
    subplot(4,3,[j j+1])
    vq = interp1(strat{i}.year + strat{i}.month/12,strat{i}.calfin_100m3,1977:1/12:2020,'linear');
    plot(1977:1/12:2020, vq)
    hold on
    h(1) = plot(strat{i}.year + strat{i}.month/12, strat{i}.calfin_100m3,'.','color',"#0072BD")
    
    vq = interp1(strat{i}.year + strat{i}.month/12,strat{i}.pseudo_100m3,1977:1/12:2020,'linear');
    plot(1977:1/12:2020, vq,'color',"#D95319")
    h(2) = plot(strat{i}.year + strat{i}.month/12, strat{i}.pseudo_100m3,'.','color',"#D95319")
    
    vq = interp1(strat{i}.year + strat{i}.month/12,strat{i}.ctyp_100m3,1977:1/12:2020,'linear');
    plot(1977:1/12:2020, vq,'color',"#EDB120")
    h(3) = plot(strat{i}.year + strat{i}.month/12, strat{i}.ctyp_100m3,'.','color',"#EDB120")

    xlim([1977 2020])
    title({titles{i},sprintf('VR = %f',VRs(i,1))})
    i=i+1;
    legend(h,'Calfin','Pseudo','Ctyp')
    ylabel('log_{10}(N+1) m^{-3}')
end
xlabel('Year')

subplot(4,3,[9 12])
plot(1:length(titles),VRs(:,1),'LineWidth',2,'Color','#a87d60')
hold on
for i=1:length(VRs)
    if VRs(i,2) == 1
        h(1)=plot(i,VRs(i,1),'*','Color','#a87d60','Linewidth',2)
    end   
end
set(gca, 'Xtick',1:length(titles),'XTickLabel',titles);
ylabel('Variance Ratio')
title('Variance Ratios')
%xlim([1 length(titles)])
legend(h,'Statistically Significant (95% CI)','Location','NW')

%% Zooplankton Transect Data

clc; clear; close all
zooplTrans = readtable('Zoo_transect_F2020.csv'); % This data is in ind m^-2

locs = ["MVCO","L01","L02","L03","L04","L05","L06","L07","L08","L09","L10","L11"];
VRs = zeros(length(locs),2);
sigTest = zeros(length(locs),2);

for i=1:length(locs)
    stat_idxs = zooplTrans(zooplTrans.Station == locs(i),:);
    stationComm = [log10(stat_idxs.CalanusFinmarchicus + 1) log10(stat_idxs.PseudocalanusMinutus+1) log10(stat_idxs.CentropagesTypicus+1)];
    VRs(i,1) = variance_ratio(stationComm);
    [lb,ub] = VR_sigTest(stationComm);
    sigTest(i,1) = lb;
    sigTest(i,2) = ub;
    if VRs(i,1) > ub || VRs(i,1) < lb
        VRs(i,2) = 1;
    end
end


figure(3)
clf
set(gcf,'color','w');
plot(1:length(locs),VRs(:,1),'LineWidth',2,'Color','#a87d60')
hold on
for i=1:length(locs)
    if VRs(i,2) == 1
        h(1)=plot(i,VRs(i,1),'*','Color','#a87d60','Linewidth',2)
    end   
end
set(gca, 'Xtick',1:length(locs),'XTickLabel',locs);
xlabel('Station')
ylabel('Variance Ratio')
title('Zooplankton Transect VR')
legend(h,'Statistically Significant (95% CI)')
xlim([1 length(locs)])



%%

function VR = variance_ratio(yy)
    % yy should be in the form of a matrix with different timeseries in each column
    % For example, columns of yy may contain the same species at different locations, or
    % different species at the same locations

    varsum = sum(nanvar(yy),'omitnan');
    covsum = sum(triu(nancov(yy),1),'all');
    VR=(varsum+2*covsum)/varsum;
end


% Significance test function
function [lb_95,ub_95] = VR_sigTest(yy) %[lb_95, ub_95]
    % Adapted from Ji's "Ji_idealized.m" reshuffling method
    % Still need to check for spurious compensatory significance results as in Solow & Dupisea (2007) 
    
    VR = variance_ratio(yy);

    [m,n]=size(yy);
    vrs_all=zeros(1000,1);
    yys=zeros(m,n);

    for i=1:1000     %permute 1000 times
        for nn=1:n           
            %random permute
            mp=randperm(m);
            yys(:,nn)=yy(mp,nn);            
        end
        
        % Calculating VR for randomly reshuffled yy timeseries       
        vrs_all(i)=variance_ratio(yys);    
    end
    
    pVal = sum(vrs_all>VR)/1000;

    mu=nanmean(vrs_all);
    ss=nanstd(vrs_all);

    % for 95% confidence bound    
    lb_95=mu-1.96*ss;   %lower bound
    ub_95=mu+1.96*ss;   %uppder bound
    
    % for 90% confidence bound    
    lb_90=mu-1.64*ss;   %lower bound
    ub_90=mu+1.64*ss;   %uppder bound
    
    % VR must be outside this limit to be significant
    
end


