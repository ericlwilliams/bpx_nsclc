function BPx_AtlasProbability

tic;
% prepare
if ~isunix,
    analy_loc = 'Z:/elw/MATLAB/bpx_analy/';
else
    analy_loc = '/Users/elw/Dropbox/eMacs/mskcc/analy/bpx_analy/';
end

fp = strcat(analy_loc,'meta_data/');
fig_loc = strcat(analy_loc,'slides/figures/latest/');

fig_ctr = 1;
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

   
do_print = true;

fn = 'BPx_DiVj_DVHs_fx-1_a2bInf.mat';dose_calib='phys';a2b_str='\infty';

fig_basename = [fig_loc, 'bpx_phys'];

screen_size=get(0,'ScreenSize');

% load data
load(strcat(fp,fn),'CGobj_current');



        LymanN = log10(CGobj_current.mLymanN);
        CGobj_current.mLymanN = LymanN;
        f = [CGobj_current.mGrp.mFlgCensor];
        
        grp = CGobj_current.mGrp;
        comps = ~f;
        eud = [grp(:).mEUD];

    CGobj_current = CGobj_current.fCrudeAtlas_relVol_DVH(-1);
    CGobj_current = CGobj_current.fBetaCumulativeProbability_DVH();

  cur_fig=figure(1); clf reset;
  %set(gcf,'Position',ss_four2three);
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGobj_current.fProbabilityFig_relVol_DVH('');
  hold on;
  
  
  set(gca,'FontSize',16)
  xlabel('Physical Dose [Gy]','FontSize',22);
  ylabel('Volume [%]','FontSize',22);
  
%  ylim([0.1,15]);
    %set(gca,'YTick',[1:2:21]); set(gca,'YTickLabel',0:10:100);
%   title('Probability of true Esophageal Toxicity $\geq$ 20\%','interpreter','latex');
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_basename,'_allfx_prob'],'-pdf');
   disp(['Saving ',fig_basename,'_allfx_prob.pdf...']);
  end

  % 5 fraction
    CGobj_current = CGobj_current.fCrudeAtlas_relVol_DVH(5);
    CGobj_current = CGobj_current.fBetaCumulativeProbability_DVH();

  cur_fig=figure(2); clf reset;
  %set(gcf,'Position',ss_four2three);
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGobj_current.fProbabilityFig_relVol_DVH('');
  hold on;
  % bin locations determined by hand
  % eg [~,med_bin] = min(abs(CGobj_current.mBinsDose - 11.7))
  h_fit=plot(83,42,'xk','MarkerSize',20,'LineWidth',4);
  h_med=plot(52,42,'+k','MarkerSize',20,'LineWidth',4);
    hold off;
  lgnd=legend([h_med h_fit],...
      ['Median split* of D$_{5\%}$ [Gy$_{3}$] for complete dataset',10,'in physical dose'],...
      'Split* at 20\% rate from regression using D$_{5\%}$ [Gy$_{3}$]',...
      'Location','North');
  set(lgnd,'FontSize',16);
  set(lgnd,'interpreter','latex');
  %ylim([0.1,15]);
  %text(0,-0.11,['*Back calculated to physical',10', 'dose delivered in 5 fractions'],...
  text(0.45,0.73,['*Back calculated to physical dose',10,'delivered in 5 fractions'],...
      'Units','normalized','interpreter','latex','fontsize',14);
   set(gca,'FontSize',16)
  xlabel('Physical Dose [Gy]','FontSize',22);
  ylabel('Volume [%]','FontSize',22);
%   title('Probability of true Esophageal Toxicity $\geq$ 20\%','interpreter','latex')
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_basename,'_5fx_prob'],'-pdf');
   disp(['Saving ',fig_basename,'_5fx_prob.pdf...']);
  end
  
  % 4 fraction
    CGobj_current = CGobj_current.fCrudeAtlas_relVol_DVH(4);
    CGobj_current = CGobj_current.fBetaCumulativeProbability_DVH();

  cur_fig=figure(3); clf reset;
  %set(gcf,'Position',ss_four2three);
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGobj_current.fProbabilityFig_relVol_DVH('');
  hold on;
  % bin locations determined by hand
  h_fit=plot(79,42,'xk','MarkerSize',20,'LineWidth',4);
  h_med=plot(50,42,'+k','MarkerSize',20,'LineWidth',4);
  hold off;
  lgnd=legend([h_med h_fit],...
      ['Median split* of D$_{5\rm{cc}}$ [Gy$_{3}$] for complete dataset',10,'in physical dose'],...
      'Split* at 20\% rate from regression using D$_{5\rm{cc}}$ [Gy$_{3}$]',...
      'Location','North');
  set(lgnd,'FontSize',16);
  set(lgnd,'interpreter','latex');
  %ylim([0.1,15]);
%   text(0,-0.11,['*Back calculated to physical',10', 'dose delivered in 4 fractions'],...
%       'Units','normalized','interpreter','latex');
 text(0.45,0.73,['*Back calculated to physical dose',10,'delivered in 4 fractions'],...
      'Units','normalized','interpreter','latex','fontsize',14);
     set(gca,'FontSize',16)
  xlabel('Physical Dose [Gy]','FontSize',22);
  ylabel('Volume [%]','FontSize',22);
%   title('Probability of true Esophageal Toxicity $\geq$ 20\%','interpreter','latex');
  
%  ylim([0.1,15]);
  
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_basename,'_4fx_prob'],'-pdf');
   disp(['Saving ',fig_basename,'_4fx_prob.pdf...']);
  end
  
    % 3 fraction
    CGobj_current = CGobj_current.fCrudeAtlas_relVol_DVH(3);
    CGobj_current = CGobj_current.fBetaCumulativeProbability_DVH();

  cur_fig=figure(4); clf reset;
  %set(gcf,'Position',ss_four2three);
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGobj_current.fProbabilityFig_relVol_DVH('');
  hold on;
  % bin locations determined by hand
  h_fit=plot(74,42,'xk','MarkerSize',20,'LineWidth',4);
   h_med=plot(47,42,'+k','MarkerSize',20,'LineWidth',4);
  hold off;
    lgnd=legend([h_med h_fit],...
      ['Median split* of D$_{5\rm{cc}}$ [Gy$_{3}$] for complete dataset',10,'in physical dose'],...
      'Split* at 20\% rate from regression using D$_{5\rm{cc}}$ [Gy$_{3}$]',...
      'Location','North');
  set(lgnd,'FontSize',16);
  set(lgnd,'interpreter','latex');

  %ylim([0.1,15]);
%   text(0,-0.11,['*Back calculated to physical',10', 'dose delivered in 3 fractions'],...
%       'Units','normalized','interpreter','latex');
   text(0.45,0.73,['*Back calculated to physical dose',10,'delivered in 3 fractions'],...
      'Units','normalized','interpreter','latex','fontsize',14);
     set(gca,'FontSize',16)
  xlabel('Physical Dose [Gy]','FontSize',22);
  ylabel('Volume [%]','FontSize',22);
%   title('Probability of true Esophageal Toxicity $\geq$ 20\%','interpreter','latex');
  
%  ylim([0.1,15]);
  
  
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_basename,'_3fx_prob'],'-pdf');
   disp(['Saving ',fig_basename,'_3fx_prob.pdf...']);
  end
  
  % 1 fraction
    CGobj_current = CGobj_current.fCrudeAtlas_relVol_DVH(1);
    CGobj_current = CGobj_current.fBetaCumulativeProbability_DVH();

  cur_fig=figure(5); clf reset;
  %set(gcf,'Position',ss_four2three);
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGobj_current.fProbabilityFig_relVol_DVH('');
  hold on;
  % bin locations determined by hand
  h_fit=plot(74,42,'xk','MarkerSize',20,'LineWidth',4);
   h_med=plot(47,42,'+k','MarkerSize',20,'LineWidth',4);
  hold off;
    lgnd=legend([h_med h_fit],...
      ['Median split* of D$_{5\rm{cc}}$ [Gy$_{3}$] for complete dataset',10,'in physical dose'],...
      'Split* at 20\% rate from regression using D$_{5\rm{cc}}$ [Gy$_{3}$]',...
      'Location','North');
  set(lgnd,'FontSize',16);
  set(lgnd,'interpreter','latex');

  %ylim([0.1,15]);
%   text(0,-0.11,['*Back calculated to physical',10', 'dose delivered in 3 fractions'],...
%       'Units','normalized','interpreter','latex');
   text(0.45,0.73,['*Back calculated to physical dose',10,'delivered in 1 fraction'],...
      'Units','normalized','interpreter','latex','fontsize',14);
     set(gca,'FontSize',16)
  xlabel('Physical Dose [Gy]','FontSize',22);
  ylabel('Volume [%]','FontSize',22);
%   title('Probability of true Esophageal Toxicity $\geq$ 20\%','interpreter','latex');
  
%  ylim([0.1,15]);
  
  
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_basename,'_1fx_prob'],'-pdf');
   disp(['Saving ',fig_basename,'_1fx_prob.pdf...']);
  end
  
    
end
    