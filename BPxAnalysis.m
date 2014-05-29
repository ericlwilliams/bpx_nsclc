function BPxAnalysis
tic;

save_result = true;

% parameters
dvhdef='DVHs';
fxnum={-1};
%fxnum={3 4 5}; % -1 all fractions, [n1 n2 ...] analyze patients with fractions of n1, n2, ... treatment planning only

%% USC BED, see C. Park et. al
% a2b_corr = 'USCBED'; 
% usc_alpha= 0.33; % Gy-1
% usc_dzero = 1.25; % Gy
% usc_dq = 1.8; % Gy
% usc_beta = ((1-usc_alpha*usc_dzero)^2)/(4*usc_dzero*usc_dq);
% usc_dtrans = (2*usc_dq)/(1-usc_alpha*usc_dzero);
% beta2alpha=[usc_beta/usc_alpha];

%% LQ Model
beta2alpha=[1/3];a2b_corr = 'BED';
%beta2alpha=[0];a2b_corr = 'PHYS';

dosestep=1;
volstep=1; % volume step in cc
timestep=3; % time step in month

% patient info.
% load patient info stored in the spread sheet

pathname='Z:/elw/MATLAB/original_data/BPX/BP_MasterDatasheet_03_07_13';

   
[~,~,xlsraw]=xlsread(pathname,'All patients');
save(pathname,'xlsraw');



% pick up related data
PtInfo = classDataFromXls();
PtInfo.mXlsRaw = xlsraw;
% MRN
%PtInfo.mLabel = 'Patient Last Name';
PtInfo.mLabel = 'MRN';
PtInfo = PtInfo.fExtractColData();
flgptcode = PtInfo.mFlg;
ptcode = PtInfo.mData;
% gender
PtInfo.mLabel = 'Sex (0-Female, 1-Male)';
PtInfo = PtInfo.fExtractColData();
flggender = PtInfo.mFlg;
ptgender = PtInfo.mData;
ptgender = strcmp(ptgender,'Male');

% Number of Fractions
PtInfo.mLabel = 'Number of Fractions';
PtInfo = PtInfo.fExtractColData();
flgfx = PtInfo.mFlg;
fx = zeros(size(flgfx));
fx(flgfx) = cell2mat(PtInfo.mData(flgfx));

% Distance to CW
PtInfo.mLabel = 'Distance to spine (cm)';
PtInfo = PtInfo.fExtractColData();
flgcm2cw = PtInfo.mFlg;
cm2cw = zeros(size(flgcm2cw));
cm2cw(flgcm2cw) = cell2mat(PtInfo.mData(flgcm2cw));

% Delivered dose
PtInfo.mLabel = 'Total Dose (cGy)';
PtInfo = PtInfo.fExtractColData();
flgtx = PtInfo.mFlg;
tx = zeros(size(flgtx));
tx(flgtx) = cell2mat(PtInfo.mData(flgtx));

% IGRT End Date
PtInfo.mLabel = 'IGRT End Date';
PtInfo = PtInfo.fExtractColData();
flgdateIGRT = PtInfo.mFlg;
dateIGRT = zeros(size(flgdateIGRT));
f = PtInfo.mData(flgdateIGRT);
dateIGRT(flgdateIGRT) = datenum(f);

% IGRT Start Date
PtInfo.mLabel = 'IGRT Start Date';
PtInfo = PtInfo.fExtractColData();
flgdateStartIGRT = PtInfo.mFlg;
dateStartIGRT = zeros(size(flgdateStartIGRT));
f = PtInfo.mData(flgdateStartIGRT);
dateStartIGRT(flgdateStartIGRT) = datenum(f);

% Date of Birth
PtInfo.mLabel = 'Date of Birth';
PtInfo = PtInfo.fExtractColData();
flgdatebirth = PtInfo.mFlg;
datebirth = zeros(size(flgdatebirth));
f = PtInfo.mData(flgdatebirth);
for i=1:length(f)
    f(i) = strrep(f(i),' ','');
end
datebirth(flgdatebirth) = datenum(f);

% age (at start of treatment)
flgAge = (flgdateStartIGRT & flgdatebirth);
age = (dateStartIGRT-datebirth)./365;

PtInfo.mLabel = 'Date of last follow-up';
PtInfo = PtInfo.fExtractColData();
flgdatefu = PtInfo.mFlg;
datefu = zeros(size(flgdatefu));
datefu(flgdatefu) = datenum(PtInfo.mData(flgdatefu));

% BPx
PtInfo.mLabel = 'Date of SBRT induced BP';
PtInfo = PtInfo.fExtractColData();
datepain = inf(size(PtInfo.mData));
f = PtInfo.mFlg;
datepain(f) = datenum(PtInfo.mData(f));

% BPx
PtInfo.mLabel = 'SBRT Induced BP (0-No, 1-Yes)';
PtInfo = PtInfo.fExtractColData();
f = PtInfo.mFlg;
censordata = PtInfo.mData(f);
censordata = cell2mat(censordata);
flgcensor = true(size(censordata));
flgcensor(censordata==1) = false;

% BPx
PtInfo.mLabel = 'Date of Recurrence Related BP';
PtInfo = PtInfo.fExtractColData();
datepain2 = inf(size(PtInfo.mData));
f = PtInfo.mFlg;
datepain2(f) = datenum(PtInfo.mData(f));

% BPx
PtInfo.mLabel = 'Reccurence Induced BP';
PtInfo = PtInfo.fExtractColData();
f = PtInfo.mFlg;
censordata = PtInfo.mData(f);
censordata = cell2mat(censordata);
flgcensor2 = true(size(censordata));
flgcensor2(censordata==1) = false;


flg = flgptcode & flgfx  & flgcm2cw  & flgdateIGRT &...
    flgdatefu & flgdatebirth & flggender &...
      flgAge & flgdateStartIGRT & flgtx;

ptcode=ptcode(flg);
% convert all mrn data to strings
containsNumbers = cellfun(@isnumeric,ptcode);
ptcode(containsNumbers) = cellfun(@num2str,ptcode(containsNumbers),...
    'UniformOutput',false);

fx=fx(flg);
tx=tx(flg);
cm2cw=cm2cw(flg);
datebirth=datebirth(flg);
dateIGRT=dateIGRT(flg);
dateStartIGRT=dateStartIGRT(flg);
datefu=datefu(flg);
datepain=datepain(flg);
datepain2=datepain2(flg);
%flgcensor=flgcensor(flg);
ptgender = ptgender(flg);
age = age(flg);




fn=['Z:\elw\MATLAB\bpx_analy\meta_data\BPx_MASTER_',dvhdef];
load(fn,'DVH');

% match patients DVH and .xls
[f1,g1]=ismember(ptcode,DVH(:,4));
[f2,g2]=ismember(DVH(:,4),ptcode);
if length(unique(g2(g2~=0)))~=length(find(g2))% || length(unique(g2~=0))~=length(find(g2))
    error('The patient ids are not unique');
end
if ~( all(f1) && all(f2) )
    disp('The patients in the spread sheet do not match with those in DVH curves, take the common part');
    % use common part only
    DVH=DVH(f2,:);
    ptcode=ptcode(f1); fx=fx(f1); tx=tx(f1);
    datebirth=datebirth(f1); cm2cw=cm2cw(f1);
    dateStartIGRT=dateStartIGRT(f1);dateIGRT=dateIGRT(f1);
    datefu=datefu(f1); 
    datepain=datepain(f1); 
    datepain2=datepain2(f1); 
    flgcensor=flgcensor(f1);
    flgcensor2=flgcensor2(f1);
    [~,g1]=ismember(ptcode,DVH(:,1));
end

CIobjs = classOutcomeIndividual();

CIobjs = classOutcomeIndividual();
CIobjs = repmat(CIobjs,[size(DVH,1),1]);
for n = 1:size(DVH,1)
    CIobjs(n,1).mID=ptcode{n};
    CIobjs(n,1).mGender = ptgender(n);
    CIobjs(n,1).mAgeAtTx = age(n);
    CIobjs(n,1).mFxNum=fx(n);
    CIobjs(n,1).mDistanceToSpine=cm2cw(n);
    CIobjs(n,1).mDoseTx=tx(n);
    CIobjs(n,1).mDosePerFx=tx(n)/fx(n);
    CIobjs(n,1).mDateBirth = datebirth(n);
    CIobjs(n,1).mDateBaseline = dateIGRT(n);
    CIobjs(n,1).mDateStartTx = dateStartIGRT(n);
    CIobjs(n,1).mDateComp = datepain(n);
    CIobjs(n,1).mDateComp2 = datepain2(n);
    CIobjs(n,1).mDateLastFollowup = datefu(n);
    CIobjs(n,1).mFlgCensor = flgcensor(n);
    CIobjs(n,1).mFlgCensor2 = flgcensor2(n);
    CIobjs(n,1).mDoseBins_org=DVH{g1(n),2}(:,1);
    CIobjs(n,1).mVolDiff=DVH{g1(n),2}(:,2);
    CIobjs(n,1).mVolCum=DVH{g1(n),2}(:,3);
    CIobjs(n,1).mBeta2AlphaCorrection = a2b_corr;    
    
     if isequal(a2b_corr,'USCBED')
        CIobjs(n,1).mAlphaUSC = usc_alpha;
        CIobjs(n,1).mDzeroUSC = usc_dzero;
        CIobjs(n,1).mDqUSC = usc_dzero;
        CIobjs(n,1).mDoseTrans = usc_dtrans;
     end
        
end

if any([CIobjs.mFxNum]==0)
    error('Fraction number is zeros');
end
CGobj_org = classOutcomeAnalysis();
CGobj_org.mLymanN = 10.^(-1:0.1:1)';

CGobj_org.mStepDose = dosestep;
CGobj_org.mStepVol = volstep;
CGobj_org.mStepTime = timestep;

CGobj_org = CGobj_org.fAddPatient(CIobjs);




% per fraction category
for n=1:length(fxnum)
    % pick up patient with the wanted fraction numbers
    if any(fxnum{n}==-1)
        f1=true(CGobj_org.mNumInGrp,1);
    else % retrieve patients with exact fractions only
        f1=arrayfun(@(x) any(x==fxnum{n}),[CGobj_org.mGrp.mFxNum]);
    end
    
    CGobj_current = CGobj_org;
    CGobj_current.mGrp = CGobj_current.mGrp(f1);
    CGobj_current.mNumInGrp = size(CGobj_current.mGrp,1);
    
   
    
    % for each alpha/beta compute a series of complication table
    disp(['Running KM and CPH analyses...']);
    for u = 1:length(beta2alpha)
        % adjust beta to alpha ratio
        
       
        
        CGobj_current.mBeta2Alpha = beta2alpha(u);
        CGobj_current = CGobj_current.fCalculateEUD();
     
        disp(['LogRankVDx_DVH...']); % V_{x} with empty patients set to zero
        CGobj_current = CGobj_current.fLogRankVDx_DVH();
        

          disp(['OveralCompCurve...']);
        CGobj_current = CGobj_current.fOverallCompCurve();
        
        disp(['Logistic Analysis']);
        CGobj_current = CGobj_current.fLogisticRegressionExact_EUD();
        
%   
        CGstrct = ObjToStruct(CGobj_current);
  
        if isequal(a2b_corr,'USCBED')
            fn=['Z:/elw/MATLAB/bpx_analy/meta_data/BPx_DiVj_',dvhdef,'_fx',num2str(fxnum{n}(1)),'_a2busc',strrep(num2str((1/beta2alpha(u)),2),'.','p'),'.mat'];
        else
         fn=['Z:/elw/MATLAB/bpx_analy/meta_data/BPx_DiVj_',dvhdef,'_fx',num2str(fxnum{n}(1)),'_a2b',num2str(1/beta2alpha(u)),'.mat'];
        end
      
        if save_result
            save(fn,'CGobj_current','CGstrct');
            disp(fn);
        end
    end
end

toc;

end