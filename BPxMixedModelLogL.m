function BPxMixedModelLogL
tic;
% prepare
fp = 'Z:\elw\MATLAB\bpx_analy\meta_data\';

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];


do_print = true;
fig_loc = 'Z:\elw\MATLAB\bpx_analy\slides\figures\latest\';

fn = {'BPx_DiVj_DVHs_fx-1_a2bInf.mat'};dose_calib='phys';
%fn = {'BPx_DiVj_DVHs_fx-1_a2b3.mat'};dose_calib='a2b3';


CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');


% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end

for m = 1:length(fn)
    CG = CGobj{m};  
    LymanN = log10(CG.mLymanN);
    CG.mLymanN = LymanN;
    
    
    % survival/complication time
    f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
    f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
    compdate = inf(CG.mNumInGrp,1);
    lastfollowup = inf(CG.mNumInGrp,1);
    compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
    lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
    compdate = min( lastfollowup, compdate );
    flgcensor = [CG.mGrp.mFlgCensor]';
    
    % Complication curve
    f1=figure(1);clf reset;hold on;
    set(f1,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    
    sa = CG.mKaplanMeierCompOverall;
    
    comp_times = sa.mSurvivalTime{1};
    prob_no_obs_sorted = sa.mSurvivalCurve{1};
    prob_obs_sorted = 1-prob_no_obs_sorted;
    norm_prob_obs_sorted = prob_obs_sorted./max(prob_obs_sorted);
    prob_no_obs_sorted = 1-norm_prob_obs_sorted;
    
    comp_times_sorted = sa.mSurvivalTimeSorted{1};
    
    prob_no_obs = -inf(size(comp_times));
    for i=1:length(comp_times)
        cur_time = comp_times(i);
        cur_times = find(cur_time<comp_times_sorted);
        if ~isempty(cur_times)
            cur_time_idx = cur_times(1);
            prob_no_obs(i) = prob_no_obs_sorted(cur_time_idx);
        else
            prob_no_obs(i) = 0;
        end
    end
    
        
    stairs(comp_times_sorted,1-prob_no_obs_sorted);
    plot(comp_times_sorted(sa.mCensorStatistics{1}(:,1)),...
        1-prob_no_obs_sorted(sa.mCensorStatistics{1}(:,1)),'+');
    
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'FontSize',12);
    xlabel('Months','fontsize',18);
    ylabel('Probability of grade >= 2 Chestwall Pain','fontsize',18);
    
    
    cur_fig=figure(2);clf reset;hold on;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    plot(comp_times,max([CG.mGrp.mDoseBins_LQ]),'*')
    set(gca,'FontSize',14);
    xlabel('Follow-up time [days]','FontSize',16);
    ylabel('D_{max} [Gy]','FontSize',16)
    
    if isequal(dose_calib,'a2b3')
            title(['BED, alpha/beta = 3 Gy'],'FontSize',18);
    else
        title(['PHYS, \alpha/\beta = \infty Gy'],'FontSize',18);
    end
        
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_dose_vs_fu_',...
            dose_calib],'-png');
        disp(['Saving ',fig_loc,'bpx_dose_vs_fu_',...
            dose_calib,'.png...']);
    end;

    
    log10a = -CG.mLymanN;
    % gEUD log-likelihood
    % using exact EUD
    euds = [CG.mGrp.mEUD]';
    pttotal = ones(CG.mNumInGrp,1);
    ptcomp = ones(CG.mNumInGrp,1); ptcomp([CG.mGrp.mFlgCensor])=0;
    loglikelihood = -inf(length(CG.mLymanN),1);
    pvals = -inf(length(CG.mLymanN),1);
    mm_loglikelihood = -inf(length(CG.mLymanN),1);
    for k=1:size(CG.mLymanN,1) % loop over each logn value, calculate fit of gEUDs
        doses=euds(:,k);
        % regression using exact EUD
        [b,~,s]=glmfit(doses,[ptcomp pttotal],'binomial','link','logit');
        pvals(k) = s.p(2);
        
        B0 = b(1);
        B1 = b(2);
        pr = exp(B0+B1*euds(:,k));
        
        pr = pr./(1+pr); % logistic probability
        pr(flgcensor) = 1-pr(flgcensor); % non-complication patients
        pr = log(pr); % log likelihood of each patients
        loglikelihood(k) = sum(pr); % loglikelihood of all
        
        %% mixed-model
        mm_pr = exp(B0+B1*euds(:,k));
        mm_pr = mm_pr./(1+mm_pr); % logistic probability
        %mm_pr(flgcensor) = ((1-mm_pr(flgcensor)).*(mm_pr(flgcensor).*prob_no_obs(flgcensor)));
        mm_pr(flgcensor) = ((1-mm_pr(flgcensor)) + (mm_pr(flgcensor).*prob_no_obs(flgcensor)));
        mm_pr = log(mm_pr);
        mm_loglikelihood(k) = sum(mm_pr);
        
        %Response functions
        cur_loga=log10a(end-k+1);
        cur_fig=figure(k+100);clf reset;%colormap(cm2);
        set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
        hold off;
        [loga,pval] = CG.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'b',1);
        CG.fComplicationObservedFig_EUD(loga,4,'b*',1);
        if loga~=cur_loga,
            disp(['loga ~= cur_loga!']);
        end
        loga_str = ['BPx',10,'Log_1_0(a) = %s',10,'p-val: %s'];
        str = sprintf(loga_str,num2str(cur_loga),num2str(pval,3));
        %text(10,0.6,str,'FontSize',16,'FontUnits','normalized');
        text(.1,0.8,str,'FontSize',16,'Units','normalized');

        
        if isequal(dose_calib,'a2b3')
                title(['BED, alpha/beta = 3 Gy'],'FontSize',18);
        else
            title(['PHYS, alpha/beta = Inf Gy'],'FontSize',18);
        end
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'bpx_logreg_response_',...
                dose_calib,'-',num2str(k)],'-png');
            disp(['Saving ',fig_loc,'bpx_logreg_response_',...
                dose_calib,'-',num2str(k),'.png...']);
        end;
    
        
%         cur_fig=figure(100+k);clf reset;hold on;
%         set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
%     
%          [rpb,rplo,rphi] = glmval(b, doses,'logit',s); % the responding function values at doses
%                 % plot
%         hold on;
%         [sorted_doses,dose_idx] = sort(doses);
%         sorted_rpb = rpb(dose_idx);
%         plot(sorted_doses,sorted_rpb,'b-','LineWidth',2);
%         %plot(doses,rpb,'b','LineWidth',2); % responding function
%         plot(doses,rpb-rplo,'b','LineWidth'); % low CI curve
%         plot(doses,rpb+rphi,'b','LineWidth'); % high CI curve
%         set(gca,'xminortick','on','yminortick','on');
%         set(gca,'box','on');
%         xlabel('gEUD'); ylabel('RP probability');
%             
    end
    
    
    % Complication curve
    cur_fig=figure(3);clf reset;hold on;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    plot(log10a,loglikelihood);
    
    
    cur_fig=figure(4);clf reset;hold on;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    plot(log10a,mm_loglikelihood,'r');
    

    cur_fig=figure(5);clf reset;hold on;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    h_ll=plot(log10a,loglikelihood,'LineWidth',2);    
    h_mm_ll=plot(log10a,mm_loglikelihood,'r','LineWidth',2);
    ylim([-39 -29]);
    xlabel('log_1_0(a)','FontSize',14);
    ylabel('Log-likelihood','FontSize',14);
    lgnd=legend([h_ll h_mm_ll],'LogL',['Mixed Model',10,'LogL'],...
        'Location','Best');
    set(lgnd,'FontSize',14);
    title('Logistic Regression Log-likelihoods','FontSize',14);
    grid on;
    
    if do_print,
       set(cur_fig,'Color','w');
       export_fig(cur_fig,[fig_loc,'logreg_llhds_',dose_calib],'-png');
    end;
    
    %% logreg pvals
    cur_fig=figure(6);clf reset;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    semilogy(log10a,pvals,'LineWidth',2);hold on;
    semilogy(log10a,repmat(0.05,length(log10a)),'r--','LineWidth',2);
    ylim([0.001 1]);
    set(gca,'FontSize',12);
    xlabel('log_1_0(a)','FontSize',14);
    ylabel('Logistic Regression p-values','FontSize',14);
    grid on;
    
    if do_print,
       set(cur_fig,'Color','w');
       export_fig(cur_fig,[fig_loc,'logreg_pvals_',dose_calib],'-png');
    end;
    
    %% Max dose logistic regression response
    
    cur_fig=figure(7);clf reset;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    
    
    mx_doses = max([CG.mGrp.mDoseBins_LQ]);
    [b,dev,st]=glmfit(mx_doses,[ptcomp pttotal],'binomial','link','logit');
    
    doses = (0:max(mx_doses))';
    [rpb,rplo,rphi] = glmval(b, doses,'logit',st); % the responding function values at mx_doses
    disp(['the beta values are: ',num2str(b')]);

    % plot
    hold on;
    plot(doses,rpb,'b','LineWidth',2); % responding function
    plot(doses,rpb-rplo,'b','LineWidth',1); % low CI curve
    plot(doses,rpb+rphi,'b','LineWidth',1); % high CI curve
    

      % grouping patients
    flg=[CG.mGrp.mFlgCensor]; % censor flags of patients
    [medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] =...
        EventObserved(flg,mx_doses,4);
    prob = numcomp./numtotal;
            % plot
    errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),'b*');
    errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),'b*');
      
    ylim([0 0.7]);
    xlims = xlim;
    xmax = xlims(2);
    xlim([0 xmax]);
    
    pval = st.p;
    pval = pval(2);
    loga_str = ['BPx',10,'Log_1_0(a) = %s (D_{max})',10,'p-val: %s'];
    str = sprintf(loga_str,'\infty',num2str(pval,3));
    %text(10,0.6,str,'FontSize',16,'FontUnits','normalized');
    text(.1,0.8,str,'FontSize',16,'Units','normalized');
    
    
    if isequal(dose_calib,'a2b3')
        title(['BED, alpha/beta = 3 Gy'],'FontSize',18);
    else
        title(['PHYS, alpha/beta = Inf Gy'],'FontSize',18);
    end
    
    
    
    set(gca,'xminortick','on','yminortick','on');
    set(gca,'box','on');
    set(gca,'FontSize',14);
    xlabel('D_{max}','FontSize',16);
    ylabel('BPx rate observed','FontSize',16);
    
    if do_print,
       set(cur_fig,'Color','w');
       export_fig(cur_fig,[fig_loc,'bpx_logreg_response_',dose_calib,'-22'],'-png');
    end;
        
end





%     euds = [CG.mGrp.mEUD]';
%     flg = [CG.mGrp.mFlgCensor]';
%
%     b0 = CG.mLogisticRegressionGridBetaRange{1};
%     b1 = CG.mLogisticRegressionGridBetaRange{2};
%
%      % for each mLymanN,b0, and b1, compute the log likelihood
%             loglikelihood = -inf(length(b0),length(b1),length(CG.mLymanN));
%             for kk = 1:length(CG.mLymanN)
%                 for jj = 1:length(b1)
%                     for ii = 1:length(b0)
%                         pr = exp(b0(ii)+b1(jj)*euds(:,kk));
%                         pr = pr./(1+pr); % logistic probability
%                         pr(flg) = 1-pr(flg); % non-complication patients
%                         pr = log(pr); % log likelihood of each patients
%                         loglikelihood(ii,jj,kk) = sum(pr); % loglikelihood of all
%                     end
%                 end
%             end
%

toc;
end
