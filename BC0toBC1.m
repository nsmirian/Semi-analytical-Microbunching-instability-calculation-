%%  from BC0 to BC1%% D2+ L1+ D3+ BC1
%%%%%%%%%%%%%%%%%%%
%% D2 Drift from BC0 to Linac 1%%
%%%%%%%%%%%%%%%%%%%%%   
%%  
LD20 = 22;            % drift path length in [m]
E0= 6.6e6;
betaxd = 15;        % average horizontal betatron function in [m]
betayd = 15;        % average vertical betatron function in [m]
% rw0 = 10e-3;      % average inner radius of round vacuum chamber in [m]
switch_bane1 =1;
sigmaz1=sigmaz0/C0;
I1=C0*I0;
%% -------------------------------% IBS %-----------------------------------%
sigdr2_ibs=IBS.drift(LD20, Efl0/me, enx, betaxd,sigmaz1 );
% IBS-induced absolute energy spread in [MeV] cumulated through d2
sigEdr2_ibs = sigdr2_ibs*Efl0
%% -------------------------------% LSC %-----------------------------------%
[S_Ldr2, ZLSCDr2]=SMI.LSCDrift( LD20,betaxd, betayd, enx,eny, Efl0, I1);
% S_Ld2 is the S-matrix^2
% ZLSCD is the longitudinal impedence 

%% %%%%%%%%%%
%% LINAC 1  %%
%%%%%%%%%%%%  
L1 = 169-123;            % Linac1 path length in [m]
Eil1= Efl0;            % initial mean energy in [eV]
Efl1 = 700e6;          % final mean energy in [eV]
betaxl1= 20;        % average horizontal betatron function in [m]
betayl1 = 15;        % average vertical betatron function in [m]
switch_bane1 =1;

%% -------------------------------% IBS %-----------------------------------%
% IBS-induced rms absolute energy spread in [MeV] with Bane's approximation
% of B-M expression, cumulated through L1
if switch_bane1 == 1
     sigL1_ibs =IBS.linac(L1,  enx, betaxl1, Eil1/me, Efl1/me,sigmaz1)
elseif switch_bane1 == 0
    sigEL1_ibs = 0;
end
%% -------------------------------% LSC %---------------------------------%

[S_L1, ZLSCL1]=SMI.LSCDlinac( L1,betaxl1, betayl1, enx,eny, Eil1,Efl1, I1)

%%
%%%%%%%%%%%%%%%%%%%
%% D3 Drift from L1 to BC1%%
%%%%%%%%%%%%%%%%%%%%%
LD3 = 182-169;            % drift path length in [m]
betaxd = 15;        % average horizontal betatron function in [m]
betayd = 15;        % average vertical betatron function in [m]
switch_bane1 =1;
%% -------------------------------% IBS %-----------------------------------%
sigdr3_ibs=IBS.drift(LD3, Efl1/me, enx, betaxd,sigmaz1)
% IBS-induced absolute energy spread in [MeV] cumulated through d2
sigEdr3_ibs = sigdr3_ibs*Efl1
%% -------------------------------% LSC %-----------------------------------%
[S_Ldr3, ZLSCDr3]=SMI.LSCDrift( LD3,betaxd, betayd, enx,eny, Efl1, I1);
% S_Ld2 is the S-matrix^2
% ZLSCD is the longitudinal impedence 
%%
%%%%%%%%%%%
%% BC1   
%%%%%%%%%%%
  %%
C1 = 8;            % linear compression factor
theta1 = 0.053233;     % dipole bending angle in [rad]
Lb1 = 0.5;        % dipole length in [m]
DL1 = 9;         % drift length between outer and inner dipole of chicane, in [m]
D112 = 2;         % drift length after dipole-2
D113 = DL1;
D114 = D112;
alpha1_in = 3;      % horizontal alpha-function at entrance of chicane
alpha1_w = 0;
beta1_in = 23;      % horizontal betatron function at entrance of chicane, in [m]
beta1_w = 10;        % horizontal betatron function at waist in the second half of chicane, in [m]
%Ef1=Efl1;
R561= -2*theta1^2*((2/3)*Lb1+DL1)    
%%%%%%%%
%%          
k1 = @(lambda) k(lambda)*C0;  
k2 = @(lambda) k1(lambda)*C1;          % compressed modulation wave number
%% -------------------------------% IBS %-----------------------------------%

sigE_l1=sqrt((sigE0_tot)^2+(sigdr3_ibs^2+sigL1_ibs^2)*Efl1^2+(sigEdr2_ibs)^2);
[sigd_BC1_w, sigdBC1_w_ibs]=IBS.DisSec(sigE_l1/Efl1,beta1_in, beta1_w, theta1, enx,...
    Efl1/me,  alpha1_in, alpha1_w , DL1,Lb1,  C1,sigmaz1)
%% -------------------------------% LSC %-----------------------------------%
S_BC1 =SMI.LSC_CSRDS(DL1,Lb1,beta1_in,beta1_w, enx, eny, Efl1, I1, theta1, ...
    C1,D112, D113, D114, sigd_BC1_w, k1)

%%
%%%%%%%%%%%%%%%%
%% MBI Matrix %%
%%%%%%%%%%%%%%%%
switch_csr=1
% Transfer matrix from injection to BC1 (included)
M_BC1 = @(lambda)S_BC1(lambda)*S_Ldr2(lambda)*S_L1(lambda)*S_Ldr3(lambda)*M_BC0(lambda);

select = @(M,row,col) M(row,col);

% Gain function after BC1
Gain1 = @(lambda) abs(select(M_BC1(lambda),1,1));

% Energy modulation induced after BC1
Emod1 = @(lambda) abs(select(M_BC1(lambda),2,1));
% Bunching factors (LSC + CSR)
abs_b1 = @(lambda) abs(select(M_BC1(lambda)*bf0(lambda),1,1))*100;

% Energy Modulation 
Emod1_keV = @(lambda) abs(select(M_BC1(lambda)*bf0(lambda),2,1))*Efl1*1e-3;
%%
% Uncorrelated energy spread RMS in induced by LSC in [keV] after BC1
sigE1_LSC = 1e-3*Efl1*sqrt(2)*sqrt((e0*c/I1)*integral(@(lambda)abs(Emod1(lambda).^2./(lambda.^2)),1e-6,200e-6,'ArrayValued',true))

% Total energy spread after BC1
sigE1_tot = sqrt((C1^2*sigd_BC1_w^2 + sigdBC1_w_ibs^2)*Efl1^2*1e-6 + ...
    sigE1_LSC^2 + switch_lh*C0^2*C1^2*sigdE_lh^2*1e-6)
%%
if fp==1
    f4=figure(4);hold on
    fplot(Emod1_keV,[1e-6,200e-6],'-b','LineWidth',3,'DisplayName','Exit of BC1')
    legend()
    set(gca,'FontSize',16)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('\DeltaE_{MBI} [keV]','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    
    f5=figure(5);hold on
    fplot(abs_b1,[1e-6,200e-6],'-b','LineWidth',3, 'DisplayName','Exit of BC1')
    set(gca,'FontSize',16)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Bunching factor [%]','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})
    legend()
%    saveas(f5,['./figures/bunching_',flagLH{switch_lh+1},'_',flagIBS{switch_bane1+1},'_',flagSpreader{switch_spreader+1},'.jpg'])
   
    f6=figure(6);hold on
    fplot(Gain1,[1e-6,200e-6],'-b','LineWidth',3,'DisplayName','Exit of BC1')
    set(gca,'FontSize',16)
    xlabel('\lambda_0 [\mum]','FontSize',16)
    ylabel('Gain','FontSize',16)
    set(gca,'XTickLabel',{'50','100','150','200'})   
    legend()
end


