function BPxDVHDisplay
tic;

% prepare
fp = 'Z:\elw\MATLAB\bpx_analy\meta_data\';

fn = {'BPx_DiVj_DVHs_fx-1_a2bInf.mat'};dose_type='phys';
%fn = {'BPx_DiVj_DVHs_fx-1_a2b3.mat'};dose_type='a2b3';


do_print = true;
fig_loc = 'Z:\elw\MATLAB\bpx_analy\slides\figures\latest\';

CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrsz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end



% plot DVHs
disp('DVH curves:');
for m = 1:length(fn)
    
    
    
    
    disp(fn{m});
    f1=figure(m); clf reset; hold on; % grid on;
    set(f1,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    f = [CGobj{m}.mGrp.mFlgCensor];
    % DVHs of censored patients
    g = find(f);
    for k = 1:length(g)
        h_cens=plot(CGobj{m}.mGrp(g(k)).mDoseBins_LQ, CGobj{m}.mGrp(g(k)).mVolCum);
    end
    % DVHs of complicated patients
    g = find(~f);
    for k = 1:length(g)
        h_comp=plot(CGobj{m}.mGrp(g(k)).mDoseBins_LQ, CGobj{m}.mGrp(g(k)).mVolCum,'r','LineWidth',2);
    end
    set(gca,'FontSize',14);
    lgnd=legend([h_comp h_cens],'With comp.','W/out comp.');
    set(lgnd,'FontSize',16);
    xlabel('Dose [Gy]','FontSize',18);
    ylabel('Vol [cc]','FontSize',18);
    set(gca,'xminortick','on','yminortick','on');
    if do_print,
       set(f1,'Color','w');
       export_fig(f1,[fig_loc,'spaghetti_',dose_type],'-png');
    end
    
    
    disp('DVH groups');
    cur_f=figure(m+1); clf reset; hold on;
    set(cur_f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    set(cur_f,'Name',fn{m});
    CGobj{m}.fDVHCurvesSummary_DVH();
    set(gca,'FontSize',14);
    xlabel('Dose [Gy]','FontSize',18);ylabel('Vol [cc]','FontSize',18);
    if do_print,
       set(cur_f,'Color','w');
       export_fig(cur_f,[fig_loc,'avg_dvhs_',dose_type],'-png');
    end
    

    %% D_{Max} histogram

    %[n_low,x_low]=hist(raw_rates(low_raw_idx==1),0:0.01:0.3);
    %bar(x_low,n_low,'g');
    
    mx = max([CGobj{m}.mGrp.mDoseBins_LQ]);
    mx_f=figure(m*10); clf reset; hold on;
    set(mx_f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    set(mx_f,'Name','D_{max}');
    if isequal(dose_type,'phys')
        mx_mx=70;
        mx_step=1;
    else % BED
        mx_mx=570;
        mx_step=10;
    end
    [n_mx,x_mx]=hist(mx,0:mx_step:mx_mx);
    [n_mx_comp,x_mx_comp]=hist(mx(~f),0:mx_step:mx_mx);
    b_mx_cens=bar(x_mx,n_mx,'b');
    b_mx_comp=bar(x_mx_comp,n_mx_comp,'r');
    lgnd_mx=legend([b_mx_cens b_mx_comp],'W/out BPx','With BPx',...
        'Location','Best');
    set(lgnd_mx,'FontSize',16);    
    set(gca,'FontSize',15);
    xlim([0 mx_mx]);
    xlabel('D_{max} [Gy]','FontSize',18);
    ylabel('Incidence','FontSize',18);
    if do_print,
       set(mx_f,'Color','w');
       export_fig(mx_f,[fig_loc,'dmax_',dose_type],'-png');
    end
           
    %% Median       
    mdn = median([CGobj{m}.mGrp.mDoseBins_LQ]);
    mdn_f=figure(m*100); clf reset; hold on;
    set(mdn_f,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    set(mdn_f,'Name','D_{median}');
   
    if isequal(dose_type,'phys')
        mx_mdn=35;
        mdn_step=1;
    else % BED
        mx_mdn=160;
        mdn_step=5;
    end
    
    [n_mdn,x_mdn]=hist(mdn,0:mdn_step:mx_mdn);
    [n_mdn_comp,x_mdn_comp]=hist(mdn(~f),0:mdn_step:mx_mdn);
    b_mdn_cens=bar(x_mdn,n_mdn,'b');
    b_mdn_comp=bar(x_mdn_comp,n_mdn_comp,'r'); 
    lgnd_mdn=legend([b_mdn_cens b_mdn_comp],'W/out BPx','With BPx',...
        'Location','Best');
    set(lgnd_mdn,'FontSize',16);
    set(gca,'FontSize',15);
    xlim([0 mx_mdn]);
    xlabel('D_{median} [Gy]','FontSize',18);
    ylabel('Incidence','FontSize',18);
    if do_print,
       set(mdn_f,'Color','w');
       export_fig(mdn_f,[fig_loc,'dmdn_',dose_type],'-png');
    end
    
end
toc;
end