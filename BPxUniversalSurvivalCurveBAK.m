function BPxUniversalSurvivalCurveBAK
tic;
% prepare
fp = 'Z:\elw\MATLAB\bpx_analy\meta_data\';

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

   
do_print = true;
fig_loc = 'Z:\elw\MATLAB\bpx_analy\slides\figures\latest\';
 

%fn = {'BPx_DiVj_DVHs_fx-1_a2busc8p6.mat'};dose_calib='uscbed';a2b_str='8p6';

fn = {'BPx_DiVj_DVHs_fx-1_a2bInf.mat'};dose_calib='phys';a2b_str='\infty';
%fn = {'BPx_DiVj_DVHs_fx-1_a2b3.mat'};dose_calib='a2b3';a2b_str='3';


CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrsz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
load(strcat(fp,fn{1}),'CGobj_current');

CGobj_current.mLymanN = 10.^(-1:0.1:1)';

 %% Set dose correction to USCBED, calculate USCBEDs
 for k=1:length([CGobj_current.mGrp])
     CGobj_current.mGrp(k).mBeta2AlphaCorrection = 'BED';
     CGobj_current.mGrp(k).mLymanN = 10.^(-1:0.1:1)';
 end

 CGobj_current.mBeta2Alpha = 1/3;
 
 %% Run basic gEUD logistic analysis first
 CGobj_lq = CGobj_current.fCalculateEUD();
        % logistic regression analysis
 CGobj_lq = CGobj_lq.fLogisticRegressionExact_EUD();

 %% printing lq
 CGobj_lq.mLymanN = log10(CGobj_lq.mLymanN);
 
 fig_ctr=1;
 cur_fig=figure(fig_ctr);  clf reset;% grid on;
 set(gcf,'Position',ss_four2three);
CGobj_lq.fLogisticRegressionLikelyhoodExactFig_a_EUD('loga','rs--',2);
  
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'lq_geud_llhds'],'-pdf');
   disp(['Saving ',fig_loc,'lq_geud_llhds.pdf...']);
 end

 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
 set(gcf,'Position',ss_four2three);
 
[loga,pval]=CGobj_lq.fLogisticRegressionRespondingCurveExactFig_a_EUD('loga','b--',2); 
CGobj_lq.fComplicationObservedFig_EUD(loga,4,'k*',2);
xlabel(['gEUD_{LQ} log_{10}(a) =',num2str(loga,1)],'FontSize',22);

 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'lq_geud_resp'],'-pdf');
   disp(['Saving ',fig_loc,'lq_geud_resp.pdf...']);
end
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
    set(gcf,'Position',ss_four2three);
 CGobj_lq.fLogisticRegressionPvalueExactFig_a_EUD('rs--',2);
  
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'lq_geud_pvals'],'-pdf');
   disp(['Saving ',fig_loc,'lq_geud_pvals.pdf...']);
end
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 
 [~,lq_aic_llhd,min_lq_aic,lq_aic_log10a] =  CGobj_lq.fLogisticRegressionAicFig_a_EUD('loga',1,'rs--',2); % 1 parameter for LQ model
  
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'lq_geud_aics'],'-pdf');
   disp(['Saving ',fig_loc,'lq_geud_aics.pdf...']);
 end

 %% Max BED logistic regression
 
 maxes = max([CGobj_lq.mGrp.mDoseBins_LQ]);
 pttotal = ones(CGobj_lq.mNumInGrp,1);
 ptcomp = ones(CGobj_lq.mNumInGrp,1); ptcomp([CGobj_lq.mGrp.mFlgCensor])=0;
 [bs,dev,st]=glmfit(maxes,[ptcomp pttotal],'binomial','link','logit');
 dpf = [dev]; % deviations
 df = [st.dfe]; % degree of freedom
 dpf = dpf./df; % deviations per degree of freedom
 loglikelyhood = -0.5*dpf;
            
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 
 doses = (0:max(maxes))'; % doses (gEUD) whose RP probability will be computed
 %doses = (0:90)';
 [rpb,rplo,rphi] = glmval(bs, doses,'logit',st); % the responding function values at doses
 %disp(['the beta values are: ',num2str(bs')]);
 
 % plot
 hold on;
 plot(doses,rpb,'b--','LineWidth',2); % responding function
 plot(doses,rpb-rplo,'b--','LineWidth',1); % low CI curve
 plot(doses,rpb+rphi,'b--','LineWidth',1); % high CI curve
 
 % grouping patients
  flg=[CGobj_lq.mGrp.mFlgCensor]; % censor flags of patients
 [medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,maxes,4);
 prob = numcomp./numtotal;
            % plot
 errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),'k*','LineWidth',2);
 errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),'k*');
 ylim([0 0.4]);
 
 set(gca,'xminortick','on','yminortick','on');
 set(gca,'box','on');
set(gca,'FontSize',18);
 textbp(['LLHD: ',num2str(loglikelyhood,3)],'fontsize',18);
 xlabel('Max BED','fontsize',22); 
 ylabel('BPx probablity','fontsize',22);
  if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'lq_max_bed_resp'],'-pdf');
   disp(['Saving ',fig_loc,'lq_max_bed_resp.pdf...']);
 end
            
%set alpha/beta AND calculate all USC BEDs
%CGobj_current.mBeta2Alpha(mUscAlpha);

%mUscAlpha = 0.052; % using 0rib cage, see http://www.rooj.com/Normal%20Tissue%20Comp.htm
mUscAlpha2Beta = 3;

%testing
%mUscRangeAlphaD0 = 0.01:0.01:0.1;
%mUscRangeDq = 0.2:0.2:1;

%best range
%mUscRangeAlphaD0 = 0.21:0.01:0.25;
%mUscRangeDq = 6.2:0.2:6.8;


%extended
%mUscRangeAlphaD0 = 0.01:0.01:0.75;
%mUscRangeAlphaD0 = 0.01:0.01:0.5;
%mUscRangeDq = 0.2:0.2:7.6;



mUscRangeAlphaD0 = 0.02:0.02:0.7;
mUscRangeDq = 0.2:0.2:15;







%% Loop over D0 Dq values, calculate DT, compute USC BED, fit outcome data           
mUscBestLog10a = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscLogLikelihoods = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscAICs = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscDt = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscPvals = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscFracLQ = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));
mUscFracFullLQ = inf(length(mUscRangeAlphaD0),length(mUscRangeDq));

for i=1:length(mUscRangeAlphaD0)
    cur_alphad0 = mUscRangeAlphaD0(i);
    for j=1:length(mUscRangeDq)
        cur_dq = mUscRangeDq(j);
        
        % calc USC BEDs
        [CGobj_current,cur_frac_lq,cur_frac_full_lq] = CGobj_current.fUSCCorrection(cur_dq,cur_alphad0,mUscAlpha2Beta);
        
        % calc gEUD
        CGobj_current = CGobj_current.fCalculateEUD();
        % logistic regression analysis
        CGobj_current = CGobj_current.fLogisticRegressionExact_EUD();
        
        %% Find best log10a and LLHDs
        st = [CGobj_current.mLogisticRegressionMat];
        dpf = [st.dev]; % deviations
        st =[st.stats];
        pvalue = [st.p];
        pvalue = pvalue(2,:); % the p-value corresponding to gEUD
        [min_p,~] = min(pvalue);
        
        df = [st.dfe]; % degree of freedom
        dpf = dpf./df; % deviations per degree of freedom
        [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
        %disp(['best log10(a) of Logistic Regression of exact gEUD is: ',num2str(-CGobj_current.mLymanN(loc))]);
        loga = -log10(CGobj_current.mLymanN(loc));
        
        %disp(['the log10(a) in coefficient searching is: ',num2str(loga)]);
        
        loglikelyhood = -0.5*dpf;
        [mxllhd,~] = max(loglikelyhood); % the maximum loglikelyhood
        
        aic = -2*loglikelyhood.*df + 2*3; % 3 parameters in usc
        [minaic,~] = min(aic);
        %disp(['Max LLHD: ',num2str(mxllhd),...
        %    ' for d0: ',num2str(cur_d0),...
        %    ' dq: ',num2str(cur_dq),...
        %    ' dt: ',num2str(cur_dt)]);
        
        mUscBestLog10a(i,j) = loga;
        mUscLogLikelihoods(i,j) = mxllhd;
        mUscAICs(i,j) = minaic;
        mUscDt(i,j) = (2*cur_dq)/(1-(cur_alphad0));
        mUscPvals(i,j) = min_p;
        mUscFracLQ(i,j) = cur_frac_lq;
        mUscFracFullLQ(i,j) = cur_frac_full_lq;
        
    end
    
    % plotyy as function of mUscDt?  will have multiple entries for each
    % mUscDt?
end
[best_llhd, mx_ind] = max(mUscLogLikelihoods(:));
[alphad0_ind, dq_ind] = ind2sub(size(mUscLogLikelihoods),mx_ind);

best_dt = mUscDt(alphad0_ind,dq_ind);
best_pval = mUscPvals(alphad0_ind,dq_ind);
best_log10a = mUscBestLog10a(alphad0_ind,dq_ind);
best_fraclq = mUscFracLQ(alphad0_ind,dq_ind);
best_fracfulllq = mUscFracFullLQ(alphad0_ind,dq_ind);
best_alphad0 = mUscRangeAlphaD0(alphad0_ind);
best_dq = mUscRangeDq(dq_ind);
disp([]);
disp([]);
disp(['== Best USC model ==',10,...
    'LLHD: ',num2str(best_llhd),10,...
    'log10a: ',num2str(best_log10a),10,...
    'pval: ',num2str(best_pval),10,...
    'alpha*D0: ',num2str(best_alphad0),10,...
    'Dq: ',num2str(best_dq),10,...
    'Dt: ',num2str(best_dt),10,...
    'Bin frac LQ: ',num2str(best_fraclq),10,...
    'Full frac LQ: ',num2str(best_fracfulllq)]);

disp([]);
disp([]);

[best_aic, min_aic_ind] = min(mUscAICs(:));
[alphad0_aic_ind, dq_aic_ind] = ind2sub(size(mUscAICs),min_aic_ind);
best_aic_llhd = mUscLogLikelihoods(alphad0_aic_ind,dq_aic_ind);
best_aic_dt = mUscDt(alphad0_aic_ind,dq_aic_ind);
best_aic_pval = mUscPvals(alphad0_aic_ind,dq_aic_ind);
best_aic_log10a = mUscBestLog10a(alphad0_aic_ind,dq_aic_ind);
best_aic_alphad0 = mUscRangeAlphaD0(alphad0_aic_ind);
best_aic_dq = mUscRangeDq(dq_aic_ind);

    
disp(['== AICs == ',10,...
    'USC min AIC: ',num2str(best_aic),10,...
     '- LLHD: ',num2str(best_aic_llhd),10,...
    '- log10a: ',num2str(best_aic_log10a),10,...
    '- pval: ',num2str(best_aic_pval),10,...
    '- alpha*D0: ',num2str(best_aic_alphad0),10,...
    '- Dq: ',num2str(best_aic_dq),10,...
    '- Dt: ',num2str(best_aic_dt)]);
disp([]);
disp(['LQ min AIC: ',num2str(min_lq_aic),10,...
    '- log10a: ',num2str(lq_aic_log10a),10,...
    '- LLHD: ',num2str(lq_aic_llhd)]);
        disp([]);
            disp([]);



 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscLogLikelihoods);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'Log-likelihood/df','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_llhds'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_llhds.pdf...']);
 end

 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscFracLQ);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'Fraction LQ Bins','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_fraclq'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_fraclq.pdf...']);
 end
 
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on; set(gcf,'Position',ss_four2three);
 set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscFracFullLQ);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'Fraction Full CW LQ ','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_fracfulllq'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_fracfulllq.pdf...']);
 end
 
 fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
 set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscAICs);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'AIC','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);

 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_aics'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_aics.pdf...']);
 end
  
 
fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
  set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscBestLog10a);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'Best log10(a)','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_log10a'],'-pdf');
    disp(['Saving ',fig_loc,'usc_geud_log10a.pdf...']);
  end
toc;

fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
  set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscDt);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'DT','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);
 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_dt'],'-pdf');
    disp(['Saving ',fig_loc,'usc_geud_dt.pdf...']);
 end
  
 
fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
  set(gcf,'Position',ss_four2three);
 imagesc(mUscRangeDq,mUscRangeAlphaD0,mUscPvals);
 set(gca,'YDir','normal');
 cb=colorbar;
 ylabel(cb,'P-value','FontSize',22);
 ylabel('\alpha\cdot D_0','FontSize',22);
xlabel('D_q [Gy]','FontSize',22);
set(gca,'FontSize',18);
xlim([min(mUscRangeDq) max(mUscRangeDq)]);
ylim([min(mUscRangeAlphaD0) max(mUscRangeAlphaD0)]);

 if do_print,
    set(gcf,'Color','w');
    export_fig(gcf,[fig_loc,'usc_geud_pvals'],'-pdf');
    disp(['Saving ',fig_loc,'usc_geud_pvals.pdf...']);
 end
  
 %% get best USC response
 % calc USC BEDs
 [CGobj_current,best_frac_lq,best_frac_full_lq] = CGobj_current.fUSCCorrection(best_dq,best_alphad0,mUscAlpha2Beta);
        
    % calc gEUD
 CGobj_current = CGobj_current.fCalculateEUD();
        % logistic regression analysis
  CGobj_current = CGobj_current.fLogisticRegressionExact_EUD();
 
  fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
 set(gcf,'Position',ss_four2three);
 
 CGobj_current.mLymanN = log10(CGobj_current.mLymanN);
    [~,pval]=CGobj_current.fLogisticRegressionRespondingCurveExactFig_a_EUD('loga','b--',2); 
    CGobj_current.fComplicationObservedFig_EUD(best_log10a,4,'k*',2);
    ylim([0 0.4]);
    
    xlabel(['gEUD_{USC} log_{10}(a) =',num2str(best_log10a,1)],'FontSize',22);

    disp([]);
disp(['Best USC response',10,...
    'p = ',num2str(pval),10,...
    'log10a = ',num2str(best_log10a),10,...
    'Frac LQ Bins= ',num2str(best_frac_lq),10,...
    'Frac Full CW LQ = ',num2str(best_frac_full_lq)]);
    disp([]);
 if do_print,
   set(cur_fig,'Color','w');
   export_fig(gcf,[fig_loc,'usc_geud_resp'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_resp.pdf...']);
 end
 

 %% profile dt
 
 dt_range = [min(mUscDt):mean(mean(diff(mUscDt))):max(mUscDt(:))]';
 dt_llhds = inf(length(dt_range),1);
 dt_pvals = inf(length(dt_range),1);
 dt_fraclq = inf(length(dt_range),1);
 dt_fracfulllq = inf(length(dt_range),1); 
 for i=1:length(dt_range)
    cur_dt = dt_range(i);
    
    [~,dt_ind] = min(abs(mUscDt(:) - cur_dt));
     
     [alphad0_dt_ind, dq_dt_ind] = ind2sub(size(mUscDt),dt_ind);
      %dt_alphad0 = mUscRangeAlphaD0(alphad0_dt_ind);     
      %dt_dq = mUscRangeDq(dq_dt_ind);     
      dt_llhds(i) = mUscLogLikelihoods(alphad0_dt_ind,dq_dt_ind);
      dt_pvals(i) = mUscPvals(alphad0_dt_ind,dq_dt_ind);
      dt_fraclq(i) = mUscFracLQ(alphad0_dt_ind,dq_dt_ind);
      dt_fracfulllq(i) = mUscFracFullLQ(alphad0_dt_ind,dq_dt_ind);
 end
 
    fig_ctr=fig_ctr+1;
    cur_fig=figure(fig_ctr);  clf reset; % grid on;
    set(gcf,'Position',ss_four2three);
 
    plot(dt_range,dt_llhds,'--','LineWidth',2);
    set(gca,'FontSize',18);
    xlabel('D_{T} [Gy]','FontSize',22);
    ylabel('Log-likelihood/df','FontSize',22);

  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'usc_geud_mxllhd_vs_dt'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_mxllhd_vs_dt.pdf...']);
  end
  
    fig_ctr=fig_ctr+1;
    cur_fig=figure(fig_ctr);  clf reset; % grid on;
    set(gcf,'Position',ss_four2three);
    
    plot(mUscRangeDq,max(mUscLogLikelihoods,[],1),'--','LineWidth',2);
    
    set(gca,'FontSize',18);
    xlabel('D_{q} [Gy]','FontSize',22);
    ylabel('Max Log-likelihood/df','FontSize',22);

  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'usc_geud_mxllhd_vs_dq'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_mxllhd_vs_dq.pdf...']);
  end
  
    
    fig_ctr=fig_ctr+1;
    cur_fig=figure(fig_ctr);  clf reset; % grid on;
    set(gcf,'Position',ss_four2three);
    
    plot(mUscRangeAlphaD0,max(mUscLogLikelihoods,[],2)','--','LineWidth',2);
    
    set(gca,'FontSize',18);
    xlabel('\alpha\cdot D_{0} [Gy]','FontSize',22);
    ylabel('Max Log-likelihood/df','FontSize',22);

  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'usc_geud_mxllhd_vs_alphad0'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_mxllhd_vs_alphad0.pdf...']);
  end
  
  fig_ctr=fig_ctr+1;
 cur_fig=figure(fig_ctr);  clf reset; % grid on;
 set(gcf,'Position',ss_four2three);
 
 
    h_full_frac_lq = plot(dt_range,dt_fracfulllq,'b--','LineWidth',2);hold on;
    h_full_bin_lq = plot(dt_range,dt_fraclq,'g--','LineWidth',2);
    hold off;
    ylim([0 1]);
   xlabel('D_{T} [Gy]','FontSize',22);

   ylabel('Fraction LQ','FontSize',22);
   lgnd = legend([h_full_frac_lq, h_full_bin_lq],'Fraction Full CW LQ','Fraction Bins LQ','location','best');
   set(lgnd,'FontSize',18);
   set(gca,'FontSize',18);
   
  if do_print,
   set(cur_fig,'Color','w');
   export_fig(cur_fig,[fig_loc,'usc_geud_lq_fracs_vs_dt'],'-pdf');
   disp(['Saving ',fig_loc,'usc_geud_lq_fracs_vs_dt.pdf...']);
 end
 
toc;

end