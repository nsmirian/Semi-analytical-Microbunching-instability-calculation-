%% from BC1 to BC2%% D4+ L2+ D5+ BC2
%%%%%%%%%%%%%%%%%%%
%% D4 Drift from BC1 to Linac 2%%
%%%%%%%%%%%%%%%%%%%%% 
I2=C1*C0*I0;
%%  
LD40 = 238-202;            % drift path length in [m]
E0= 6.6e6;
betaxd = 10;        % average horizontal betatron function in [m]
betayd = 10;        % average vertical betatron function in [m]
% rw0 = 10e-3;        % average inner radius of round vacuum chamber in [m]
switch_bane1 =1;
sigmaz2=sigmaz1/C1;
I2=C1*I1;
%% -------------------------------% IBS %-----------------------------------%
sigdr4_ibs=IBS.drift(LD40, Efl1/me, enx, betaxd,sigmaz2);
% IBS-induced absolute energy spread in [MeV] cumulated through d2
sigEdr4_ibs = sigdr4_ibs*Efl1
%% -------------------------------% LSC %-----------------------------------%
[S_Ldr4, ZLSCDr4]=SMI.LSCDrift( LD40,betaxd, betayd, enx,eny, Efl1, I2);
%% %%%%%%%%%%
%% LINAC 2  %%
%%%%%%%%%%%%  
L2 = 381-238;            % Linac1 path length in [m]  73 m
Eil2= Efl0;            % initial mean energy in [eV]
Efl2 = 2500e6;          % final mean energy in [eV]
betaxl2= 20;        % average horizontal betatron function in [m]
betayl2 = 15;        % average vertical betatron function in [m]
switch_bane1 =1;

%% -------------------------------% IBS %-----------------------------------%
% IBS-induced rms absolute energy spread in [MeV] with Bane's approximation
% of B-M expression, cumulated through L1
if switch_bane1 == 1
     sigdL2_ibs =IBS.linac(L2,  enx, betaxl2, Eil2/me, Efl2/me,sigmaz2)
elseif switch_bane1 == 0
    sigEL2_ibs = 0;
end
%% -------------------------------% LSC %---------------------------------%

[S_L2, ZLSCL2]=SMI.LSCDlinac( L2,betaxl2, betayl2, enx,eny, Eil2,Efl2, I2)

%%
%%%%%%%%%%%%%%%%%%%
%% D5 Drift from L2 to BC2%%
%%%%%%%%%%%%%%%%%%%
LD5 = 393-381;            % drift path length in [m]
betaxd = 10;        % average horizontal betatron function in [m]
betayd = 10;        % average vertical betatron function in [m]
switch_bane1 =1;
%% -------------------------------% IBS %-----------------------------------%
sigdr5_ibs=IBS.drift(LD5, Efl2/me, enx, betaxd,sigmaz2)
% IBS-induced absolute energy spread in [MeV] cumulated through d2
sigEdr5_ibs = sigdr5_ibs*Efl2
%% -------------------------------% LSC %-----------------------------------%
[S_Ldr5, ZLSCDr5]=SMI.LSCDrift( LD5,betaxd, betayd, enx,eny, Efl2, I2);
% S_Ld2 is the S-matrix^2
% ZLSCD is the longitudinal impedence 
%%
%%%%%%%%%%%
% BC2   
%%%%%%%%%%%
  %%
C2 = 5;            % linear compression factor
theta2 = 0.0412;     % dipole bending angle in [rad]
Lb2 = 0.5;        % dipole length in [m]
DL2 = 402-393;         % drift length between outer and inner dipole of chicane, in [m]
D212 =404-402;         % drift length after dipole-2
D213 = DL2;
D214 = D212;
alpha2_in = 3;      % horizontal alpha-function at entrance of chicane
alpha2_w = 0;
beta2_in = 23;      % horizontal betatron function at entrance of chicane, in [m]
beta2_w = 3;        % horizontal betatron function at waist in the second half of chicane, in [m]
%Ef1=Efl1;
R562= -2*theta2^2*((2/3)*Lb2+DL2)    
%%%%%%%%
%%          
k1 = @(lambda) k(lambda)*C0;  
k2 = @(lambda) k1(lambda)*C1;          % compressed modulation wave number
k3 = @(lambda) k2(lambda)*C2;
%% -------------------------------% IBS %-----------------------------------%
sigE_l2=sqrt((sigE1_tot)^2+(sigdr4_ibs^2+sigdL2_ibs^2+sigdr5_ibs^2)*Efl2^2);
[sigd_BC2_w, sigdBC2_w_ibs]=IBS.DisSec(sigE_l2/Efl2,beta2_in, beta2_w, theta2, enx,...
    Efl2/me,  alpha2_in, alpha2_w , DL2,Lb2,  C2,sigmaz2);
%% -------------------------------% LSC %-----------------------------------%
S_BC2 =SMI.LSC_CSRDS(DL2,Lb2,beta2_in,beta2_w, enx, eny, Efl2, I2, theta2, ...
    C2,D212, D213, D214, sigd_BC2_w, k2);
%%

%%%%%%%%%%%%%%%%
%% MBI Matrix %%
%%%%%%%%%%%%%%%%
switch_csr=1
% Transfer matrix from injection to BC1 (included)
M_BC2 = @(lambda)S_BC2(lambda)*S_Ldr5(lambda)*S_L2(lambda)*S_Ldr4(lambda)*M_BC1(lambda);

select = @(M,row,col) M(row,col);

% Gain function after BC2
Gain2 = @(lambda) abs(select(M_BC2(lambda),1,1));

% Energy modulation induced after BC2
Emod2 = @(lambda) abs(select(M_BC2(lambda),2,1));
% Bunching factors (LSC + CSR)
abs_b2 = @(lambda) abs(select(M_BC2(lambda)*bf0(lambda),1,1))*100;

% Energy Modulation 
Emod2_keV = @(lambda) abs(select(M_BC2(lambda)*bf0(lambda),2,1))*Efl2*1e-3;
%%
% Uncorrelated energy spread RMS in induced by LSC in [keV] after BC2
sigE2_LSC = 1e-3*Efl2*sqrt(2)*sqrt((e0*c/I2)*integral(@(lambda)abs(Emod2(lambda).^2./(lambda.^2)),1e-4,200e-4,'ArrayValued',true))

% Total energy spread after BC2
sigE2_tot = sqrt((C2^2*sigd_BC2_w^2 + sigdBC2_w_ibs^2)*Efl2^2*1e-6 +...
    sigE2_LSC^2 + switch_lh*C0^2*C1^2*C2^2*sigdE_lh^2*1e-6)
%%
% f4=figure(4);hold on
%     fplot(Emod2_keV,[1e-6,200e-6],'-g','LineWidth',3, 'DisplayName','Exit of BC2')
%     set(gca,'FontSize',16)
%     xlabel('\lambda_0 [\mum]','FontSize',16)
%     ylabel('\DeltaE_{MBI} [keV]','FontSize',16)
%     set(gca,'XTickLabel',{'50','100','150','200'})
%      legend()
%      
% f5=figure(5);hold on
%     fplot(abs_b2,[1e-6,200e-6],'-g','LineWidth',3, 'DisplayName','Exit of BC2')
%     set(gca,'FontSize',16)
%     xlabel('\lambda_0 [\mum]','FontSize',16)
%     ylabel('Bunching factor [%]','FontSize',16)
%     set(gca,'XTickLabel',{'50','100','150','200'})
%     legend()
% %    saveas(f5,['./figures/bunching_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])
%    
%  f6=figure(6);hold on
%     fplot(Gain2,[1e-6,200e-6],'-g','LineWidth',3, 'DisplayName','Exit of BC2')
%     set(gca,'FontSize',16)
%     xlabel('\lambda_0 [\mum]','FontSize',16)
