%%%%%%%%%%%%%%%%%%%%%
%% D0  First Drift %%
%%%%%%%%%%%%%%%%%%%%%
run('constants.m')
run('functions.m')
LD = 25;            % drift path length in [m]
E0= 6.6e6;          % initial and final mean energy in [eV]
betaxd = 10;        % average horizontal betatron function in [m]
betayd = 10;        % average vertical betatron function in [m]
rw0 = 10e-3;        % average inner radius of round vacuum chamber in [m]
switch_bane1 =0;
%% IBS& LSC
%-------------------------------% IBS %-----------------------------------%
sigdr0_ibs=0; %sigddrift_ibs(LD, E0/me, enx, betaxd);Energy is so low so forget IBS
% IBS-induced absolute energy spread in [MeV] cumulated through d2
sigEr0_ibs = sigdr0_ibs*E0;
%-------------------------------% LSC %-----------------------------------%
[S_dr0, ZLSCDr0]=SMI.LSCDrift( LD,betaxd, betayd, enx,eny, E0, I0);
% S_dr0 is the S-matrix^2
% ZLSCD0 is the longitudinal impedence 
%%%%%%%%%%%%
%% LINAC0 %%
%%%%%%%%%%%%  
L0 = 20;            % Linac1 path length in [m]
Eil0= E0;            % initial mean energy in [eV]
Efl0 = 126.3e6;          % final mean energy in [eV]
betaxl0= 20;        % average horizontal betatron function in [m]
betayl0 = 15;        % average vertical betatron function in [m]
switch_bane1 =0;
sigmaz0 = c*Q/(5.5*I0);  % buncg length
% rw0 = 10e-3;       % average inner radius of round vacuum chamber in [m]
%% 
%-------------------------------% IBS %-----------------------------------%
sigd_L0 = sqrt(sigd0^2 + sigdr0_ibs^2);  % Uncorrelated energy spread at the entrance of L1
% IBS-induced rms absolute energy spread in [MeV] with Bane's approximation
if switch_bane1 == 1
     sigL0_ibs =IBS.drift (L0,  enx, betaxl0, Eil0/me, Efl0/me,sigmaz0)
elseif switch_bane1 == 0
    sigL0_ibs = 0;
end
%-------------------------------% LSC %-----------------------------------%
[S_L0, ZLSCL0]=SMI.LSCDlinac( L0,betaxl0, betayl0, enx,eny, Eil0,Efl0, I0)
%%
%%%%%%%%%%%%%%%%%%%
%%  D2=Drift +LH 
%%%%%%%%%%%%%%%%%%%
LD0 = 27;            % drift path length in [m]
% Eid0 =Efl0 ;           % initial mean energy in [eV]
% Efd0 = Efl0 ;           % final mean energy in [eV]
betaxd0 = 10;        % average horizontal betatron function in [m]
betayd0 = 10;        % average vertical betatron function in [m]
rw0 = 10e-3;        % average inner radius of round vacuum chamber in [m]
switch_bane1 =1
                                                                % compressed modulation wave number
%%
%-------------------------------% IBS %-----------------------------------%
sigdr1_ibs=IBS.drift(LD0, Efl0/me, enx, betaxd0,sigmaz0);
% IBS-induced absolute energy spread in [MeV] cumulated through d2
sigEdr1_ibs = sigdr1_ibs*Efl0;
%-------------------------------% LSC %-----------------------------------%
[S_dr1, ZLSCDr1]=SMI.LSCDrift( LD0,betaxd0, betayd0, enx,eny, Efl0, I0);
%%
sigd_dr1 = sqrt((sigd0^2 + sigdr0_ibs^2)*E0^2+ (sigdr1_ibs^2+sigL0_ibs^2)*Efl0^2)
%%
%%%%%%%%%%%
%% Dogleg =Spreader   At this script I considered dogleg as a drift :) 
%%%%%%%%%%%
L_dog=95-73;
LD =L_dog;            % drift path length in [m]
En=Efl0               % initial mean energy in [eV]
betaxd = 10;          % average horizontal betatron function in [m]
betayd = 10;          % average vertical betatron function in [m]
rw0 = 10e-3;          % average inner radius of round vacuum chamber in [m]
switch_bane1 =1;
%-------------------------------% IBS %-----------------------------------%
if switch_bane1 == 1
    sigd_dog_ibs=IBS.drift(LD, En/me, enx, betaxd,sigmaz0);
% IBS-induced absolute energy spread in [MeV] cumulated through 
elseif switch_bane1 == 0
    sigd_dog_ibs = 0;
end
sigEd2_ibs = sigd_dog_ibs*En;
%-------------------------------% LSC %-----------------------------------%
[S_L_dog, ZLSC_dog]=SMI.LSCDrift( LD,betaxd, betayd, enx,eny, En, I0); 
%%
%%
sigE_dog=sqrt((sigd0^2+sigdr0_ibs^2)*E0^2+(sigdr1_ibs^2+sigL0_ibs^2+sigd_dog_ibs^2)*Efl0^2)  %ev RMS energy spread
%%%%%%%%%%%
%% BC0   
%%%%%%%%%%%
C0 = 3;            % linear compression factor
theta0 = 0.136659;     % dipole bending angle in [rad]
Lb0 = 0.5;        % dipole length in [m]
DL0 = 1.0;         % drift length between outer and inner dipole of chicane, in [m]
D012 = 1.5;         % drift length after dipole-2
D013 = DL0;
D014 = D012;
alpha0_in = 3;      % horizontal alpha-function at entrance of chicane
alpha0_w = 0;
beta0_in = 23;      % horizontal betatron function at entrance of chicane, in [m]
beta0_w = 3;        % horizontal betatron function at waist in the second half of chicane, in [m]
%Ef0=Efl0;
R560 = -2*theta0^2*((2/3)*Lb0+DL0);    
%%%%%%%%
%%                                              
k1 = @(lambda) k(lambda)*C0;     % compressed modulation wave number
%% ّIBS &LSC & CSR--------------------------------------------------------%
%-------------------------------% IBS %-----------------------------------%
[sigd_BC0_w, sigdBC0_w_ibs]=IBS.DisSec(sigE_dog/Efl0,beta0_in, beta0_w, theta0, enx,...
    Efl0/me,  alpha0_in, alpha0_w , DL0,Lb0,  C0, sigmaz0)
%%-------------------------------% LSC %-----------------------------------%
S_BC0 =SMI.LSC_CSRDS(DL0,Lb0,beta0_in,beta0_w, enx, eny, Efl0, I0, theta0, ...
    C0,D012, D013, D014, sigd_BC0_w, k)
%%
close all
%%%%%%%%%%%%%%%%
%% MBI Matrix %%
%%%%%%%%%%%%%%%%
% Transfer matrix from injection to BC0 (included)
M_BC0 = @(lambda) S_BC0(lambda)*S_L_dog(lambda)*S_dr1(lambda)*S_L0(lambda)*S_dr0(lambda);

select = @(M,row,col) M(row,col);
% Gain function after D0+L0+BC0
Gain0 = @(lambda) abs(select(M_BC0(lambda),1,1));

% Energy modulation induced after D0+L0+BC0
Emod0 = @(lambda) abs(select(M_BC0(lambda),2,1));
% Bunching factors (LSC + CSR)
abs_b0 = @(lambda) abs(select(M_BC0(lambda)*bf0(lambda),1,1))*100;

% Energy Modulation 
Emod0_keV = @(lambda) abs(select(M_BC0(lambda)*bf0(lambda),2,1))*Efl0*1e-3;

% Uncorrelated energy spread RMS in induced by LSC in [keV] after BC0
sigE0_LSC = 1e-3*Efl0*sqrt(2)*sqrt((e0*c/I0)*integral(@(lambda)abs(Emod0(lambda).^2./(lambda.^2)),1e-6,200e-6,'ArrayValued',true));

% Total energy spread after BC0
sigE0_tot = sqrt((C0^2*sigd_BC0_w^2 + sigdBC0_w_ibs^2)*Efl0^2*1e-6 + sigE0_LSC^2 + switch_lh*C0^2*sigdE_lh^2*1e-6);

%%
if fp==1
    f4=figure(14);
    fplot(Emod0_keV,[1e-6,200e-6],'-r','LineWidth',3,'DisplayName','Exit of BC0')
    legend()
    set(gca,'FontSize',16)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('\DeltaE_{MBI} [keV]','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    
    f5=figure(15);
    fplot(abs_b0,[1e-6,200e-6],'-r','LineWidth',3, 'DisplayName','Exit of BC0')
    set(gca,'FontSize',16)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Bunching factor [%]','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    legend()
    %    saveas(f5,['./figures/bunching_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])
    
    f6=figure(16);
    fplot(Gain0,[1e-6,200e-6],'-r','LineWidth',3,'DisplayName','Exit of BC0')
    set(gca,'FontSize',16)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Gain','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    legend()
end
    
