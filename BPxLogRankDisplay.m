function BPxLogRankDisplay
tic;
% prepare
if ~isunix,
    analy_loc = 'Z:/elw/MATLAB/bpx_analy/';
else
    analy_loc = '/Users/elw/Dropbox/eMacs/mskcc/analy/bpx_analy/';
end

fp = strcat(analy_loc,'meta_data/');
fig_loc = strcat(analy_loc,'slides/figures/latest/');

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

   
do_print = true;

 

%fn = {'BPx_DiVj_DVHs_fx-1_a2busc8p6.mat'};dose_calib='uscbed';a2b_str='8p6';

%fn = {'BPx_DiVj_DVHs_fx-1_a2bInf.mat'};dose_calib='phys';a2b_str='\infty';
fn = {'BPx_DiVj_DVHs_fx-1_a2b3.mat'};dose_calib='a2b3';a2b_str='3';


CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrsz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end

for m = 1:length(fn)
    
    
    %% KM curve displays
    
    % compose complication KM curves
    
    
    
    CG = CGobj{m};
    % survival/complication time
    f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
    f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
    compdate = inf(CG.mNumInGrp,1);
    lastfollowup = inf(CG.mNumInGrp,1);
    compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
    lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
    compdate = min( lastfollowup, compdate );
    flgcensor = [CG.mGrp.mFlgCensor]';
    
    
    
    %% Compare SBRT BPx and Recurrence BPx
    f2 = ~cellfun('isempty',{CG.mGrp.mDateComp2}); % patients with no complication date
    compdate2 = inf(CG.mNumInGrp,1);
    compdate2(f2) = ([CG.mGrp(f2).mDateComp2] - [CG.mGrp(f2).mDateBaseline])' / 30;
    compdate2 = min( lastfollowup, compdate );
    flgcensor2 = [CG.mGrp.mFlgCensor2]';
    
    sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
    % assign properties of object sa
    
    survivedate={compdate; compdate2}; % survive time of each group
    fcensor={flgcensor; flgcensor2}; % censor flag for each group
    sa.mSurvivalTime=survivedate;
    sa.mFlgCensor=fcensor;
    % compute survival curves and compare them
    sa=sa.fCalculateSurvivalCurve();
    sa=sa.fCombineSurvivalTime();
    sa=sa.fCompareSurvivalByLogrank();
    
    
    
    % plot
    figure(1);  clf reset; hold on; % grid on;
    set(gcf,'Position',ss_four2three);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    
    str_pval1 = {strcat('Log-Rank p-value = ',num2str(sa.mpValue,2))};
    text(84,0.03,str_pval1,'FontSize',14);
    lgnd=legend(h_km,...
        strcat('SBRT BPx'),...
        strcat('Recurrence BPx'),...
        'Location','SouthEast');
    set(lgnd,'FontSize',16);
    %h=legend;
    %set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',14);
    ylabel(['Probability of BPx'],'fontsize',14);
    title('Brachial Plexopathy Incidence','fontsize',14);
    
    %% Prescribed dose
    tx = [CG.mGrp.mDoseTx];
    med_tx = median(tx);
    flg_below=tx<=med_tx;
    
    survivedate={compdate(flg_below); compdate(~flg_below)}; % survive time of each group
    fcensor={flgcensor(flg_below); flgcensor(~flg_below)}; % censor flag for each group
    sa.mSurvivalTime=survivedate;
    sa.mFlgCensor=fcensor;
    % compute survival curves and compare them
    sa=sa.fCalculateSurvivalCurve();
    sa=sa.fCombineSurvivalTime();
    sa=sa.fCompareSurvivalByLogrank();
    
    % plot
    figure(2);  clf reset; hold on; % grid on;
    set(gcf,'Position',ss_four2three);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    
%     str_pval1 = {strcat('Log-Rank p-value = ',num2str(sa.mpValue,2))};
    text(84,0.03,str_pval1,'FontSize',14);
    lgnd=legend(h_km,...
        ['Tx $\leq',num2str(med_tx),'$'],...
        ['Tx $> ',num2str(med_tx),'$'],...
        'Location','SouthEast');
    set(lgnd,'FontSize',16);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',14);
    ylabel(['Probability of BPx'],'fontsize',14);
    title('Brachial Plexopathy Incidence','fontsize',14);
    
      %% Max BED Dose
    mx = max([CG.mGrp.mDoseBins_LQ]);
    med_mx = median(mx);
    flg_below=mx<=med_mx;
    
    survivedate={compdate(flg_below); compdate(~flg_below)}; % survive time of each group
    fcensor={flgcensor(flg_below); flgcensor(~flg_below)}; % censor flag for each group
    sa.mSurvivalTime=survivedate;
    sa.mFlgCensor=fcensor;
    % compute survival curves and compare them
    sa=sa.fCalculateSurvivalCurve();
    sa=sa.fCombineSurvivalTime();
    sa=sa.fCompareSurvivalByLogrank();
    
    % plot
    figure(3);  clf reset; hold on; % grid on;
    set(gcf,'Position',ss_four2three);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1}./12,1-sa.mSurvivalCurve{1},'LineWidth',2);
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+','MarkerSize',15);
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2}./12,1-sa.mSurvivalCurve{2},'r','LineWidth',2);
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1))./12,...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+','MarkerSize',15);
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    
        % calc. hazard ratio from split
    cox_beta=coxphfit(~flg_below',compdate,'baseline',0,'censoring',flgcensor);
    cur_lr_hrs = exp(cox_beta);
    str_hr = ['HR = ',num2str(cur_lr_hrs,3)];

    ylim([0 0.2]);
%    set(gca,'FontSize',16);
    str_pval1 = {strcat('Log-Rank p-value = ',num2str(sa.mpValue,'%1.2g'))};
    %text(75,0.03,[str_pval1,str_hr],'FontSize',18);
    textbp([str_pval1,str_hr],'FontSize',18);
    
    lgnd=legend(h_km,...
        ['D$_{\rm{max}} \leq',num2str(med_mx,3),'$~Gy$_{',a2b_str,'}$'],...
        ['D$_{\rm{max}} >',num2str(med_mx,3),'$~Gy$_{',a2b_str,'}$'],...
        'Location','Best');
    set(lgnd,'FontSize',18);
    h=legend;
    set(h,'interpreter','latex');
    
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Years'],'fontsize',22);
    ylabel(['Probability of BPx'],'fontsize',22);
    set(gca,'FontSize',18)
    %title('Brachial Plexopathy Incidence','fontsize',20);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'bpx_km_',dose_calib,'_dmax'],'-pdf');
   disp(['Saving ',fig_loc,'bpx_km_',dose_calib,'_dmax.pdf']);
 end
%         %% BED D05
    doses = [CG.mGrp.mDoseBins_LQ];
    vols = [CG.mGrp.mVolCum];
    d05s = zeros(size(vols,2),1);
    for j=1:size(vols,2) % for each patient
        vol = vols(:,j);
        vol = vol./max(vol);
        v05 = vol<0.05;
        d05_inds = find(v05);
        dose = doses(:,j);
        d05 = min(dose(d05_inds));
        d05s(j) = d05;
    end
    
      
    med_d05s = median(d05s);
    flg_below=d05s<=med_d05s;
    
    survivedate={compdate(flg_below); compdate(~flg_below)}; % survive time of each group
    fcensor={flgcensor(flg_below); flgcensor(~flg_below)}; % censor flag for each group
    sa.mSurvivalTime=survivedate;
    sa.mFlgCensor=fcensor;
    % compute survival curves and compare them
    sa=sa.fCalculateSurvivalCurve();
    sa=sa.fCombineSurvivalTime();
    sa=sa.fCompareSurvivalByLogrank();
    
    % plot
    figure(4);  clf reset; hold on; % grid on;
    set(gcf,'Position',ss_four2three);
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1}./12,1-sa.mSurvivalCurve{1},'LineWidth',2);
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+','MarkerSize',15);
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2}./12,1-sa.mSurvivalCurve{2},'r','LineWidth',2);
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1))./12,...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+','MarkerSize',15);
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));

    % calc. hazard ratio from split
    cox_beta=coxphfit(~flg_below,compdate,'baseline',0,'censoring',flgcensor);
    cur_lr_hrs = exp(cox_beta);
    
    str_hr = ['HR = ',num2str(cur_lr_hrs,3)];
    
    ylim([0 0.2]);
    
    set(gca,'FontSize',18);

    str_pval1 = {strcat('Log-Rank p-value = ',num2str(sa.mpValue,'%3.2g'))};
    %text(75,0.03,[str_pval1,str_hr],'FontSize',18);
    textbp([str_pval1,str_hr],'FontSize',18);
    lgnd=legend(h_km,...
        ['D$_{05} \leq',num2str(med_d05s,3),'$~Gy$_{',a2b_str,'}$'],...
        ['D$_{05} >',num2str(med_d05s,3),'$~Gy$_{',a2b_str,'}$'],...
        'Location','best');
    set(lgnd,'FontSize',18);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Years'],'fontsize',22);
    ylabel(['Probability of BPx'],'fontsize',22);
    %title('Brachial Plexopathy Incidence','fontsize',20);
    
   if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'bpx_km_',dose_calib,'_d05'],'-pdf');
   disp(['Saving ',fig_loc,'bpx_km_',dose_calib,'_d05.pdf']);
 end
    %% EUDs
    euds = [CG.mGrp.mEUD];
    lyman_n = CG.mLymanN;
    log10a = log10([CG.mLymanN].^-1);
    
    lr_pvals = inf(size(log10a));
    lr_hrs = inf(size(log10a));
    lr_hrs_ucl = inf(size(log10a));
    lr_hrs_lcl = inf(size(log10a));
    
    cur_fig_ctr=4;
    for i=size(euds,1):-1:1
        cur_fig_ctr = cur_fig_ctr+i;
        eud = euds(i,:);
        med_eud = median(eud);
        flg_below=eud<=med_eud;
        
        survivedate={compdate(flg_below); compdate(~flg_below)}; % survive time of each group
        fcensor={flgcensor(flg_below); flgcensor(~flg_below)}; % censor flag for each group
        sa.mSurvivalTime=survivedate;
        sa.mFlgCensor=fcensor;
        % compute survival curves and compare them
        sa=sa.fCalculateSurvivalCurve();
        sa=sa.fCombineSurvivalTime();
        sa=sa.fCompareSurvivalByLogrank();
        lr_pvals(i)=sa.mpValue;
        
        disp([num2str(log10a(i)),' p: ',num2str(lr_pvals(i)),' areas: ',...
            num2str(sa.mCurveArea(1)),', ',num2str(sa.mCurveArea(2))]);
        % calc. hazard ratio from split
        [cox_beta,~,~,cox_stats]=coxphfit(~flg_below',compdate,'baseline',0,'censoring',flgcensor);
        lr_hrs(i) = exp(cox_beta);
        lr_hrs_ucl(i) = lr_hrs(i) + 1.96*cox_stats.se;
        lr_hrs_lcl(i) = lr_hrs(i) - 1.96*cox_stats.se;
        
        % plot
        cur_fig=figure(cur_fig_ctr);  clf reset; hold on; % grid on;
        set(gcf,'Position',ss_four2three);
        h_km(1)=stairs(sa.mSurvivalTimeSorted{1}./12,1-sa.mSurvivalCurve{1},'LineWidth',2);
        plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1))./12,...
            1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+','MarkerSize',15);
        h_km(2)=stairs(sa.mSurvivalTimeSorted{2}./12,1-sa.mSurvivalCurve{2},'r','LineWidth',2);
        plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1))./12,...
            1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+','MarkerSize',14);
        %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
        ylim([0 0.2]);
        %xlim([0 110]);
        str_pval1 = ['p-value = ',num2str(sa.mpValue,2)];
        %text(84,0.03,str_pval1,'FontSize',14);
        str_hr = ['HR = ',num2str(lr_hrs(i),3)];
        str_log10a=['log_1_1(a) = ',num2str(log10a(i))];
        text(84,0.03,[str_log10a,10,str_pval1,10,str_hr],'FontSize',16);
        lgnd=legend(h_km,...
           ['EUD <= ',num2str(med_eud,2)],...
            ['EUD > ',num2str(med_eud,2)],...
            'Location','NorthEast');
      set(lgnd,'FontSize',18);
 %       h=legend;
 %       set(h,'interpreter','latex');
  
        set(gca,'FontSize',18);
        set(gca,'xminortick','on','yminortick','on');
        xlabel(['Years'],'fontsize',22);
        ylabel(['Probability of BPx'],'fontsize',22);
        title('Brachial Plexopathy Incidence','fontsize',22);
        
       if do_print,
           %saveas(cur_fig,[fig_loc,'km_bpx_log10a-',num2str(i),'.eps']);
           set(cur_fig,'Color','w');
           export_fig(cur_fig,[fig_loc,'km_bpx_log10a_',dose_calib,'-',num2str(i)],...
               '-png');
           %print_fig(cur_fig,fig_loc,['km_bpx_log10a-',num2str(i)],'pdf');
       end;
    end
end

    cur_fig_ctr=cur_fig_ctr+1;
   f=figure(cur_fig_ctr); clf reset;
   
    set(f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    
    x_axis = log10a;
    [ax,h1,h2]=plotyy(x_axis,lr_pvals,x_axis,lr_hrs,@semilogy);hold on;
    set(h1,'LineWidth',2);
    set(h2,'LineWidth',2);
    set(ax(1),'YLim',[0.001 1]);
    set(get(ax(1),'Ylabel'),'String','p-value','FontSize',20);    

    set(ax(2),'YLim',[1 10]);
    set(get(ax(2),'Ylabel'),'String','Hazard Ratio','FontSize',20);
    semilogy(x_axis,repmat(0.05,length(x_axis)),'r--','LineWidth',1);
    set(ax(1),'XLim',[min(log10a) max(log10a)]);
    set(ax(2),'XLim',[min(log10a) max(log10a)]);
    set(ax,'FontSize',15);
 
    hold off; % grid on;
    xlabel('log_1_0(a)','fontsize',18);
    %set(lgnd,'FontSize',14);
    
 if do_print,
           %saveas(cur_fig,[fig_loc,'km_bpx_log10a-',num2str(i),'.eps']);
           set(f,'Color','w');
           export_fig(cur_fig,[fig_loc,'bpx_logrank_pval_hr_',dose_calib],...
               '-pdf');
           %print_fig(cur_fig,fig_loc,['km_bpx_log10a-',num2str(i)],'pdf');
       end;

toc;
end