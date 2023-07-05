%% Constant
A_csr = 1.63-1i*0.9454;   % CSR free vacuum impedance
Z0 = 376.73;              % vacuum impedance n [Ohm]
e0 = 1.6e-19;              % electron charge in [C]
me = 0.511e6;               % electron rest mass in [eV]
re = 2.818e-15;           % electron classical radius in [m]
IA = 17045;               % Alfven current in [A]
c = 2.998e8;              % light speed in vacuum in [m/s]


%%
%%%%%%%%%%%%%%
% SWITCHES  %%
%%%%%%%%%%%%%%

switch_lh = 1;            % 1 means heating depends on transverse size, 0 means heating is purely Gaussian
switch_csr = 1;           % 1 means CSR ON, 0 means CSR OFF
switch_cer = 1;           % 1 means CER ON, 0 means CER OFF
switch_bane1 = 1;         % 1 means IBS effects ON, 0 means IBS effects OFF
switch_FEL = 1;           % 1 menas FEL2, 0 means FEL1
switch_spreader = 1;      % 3 means Wolski matrix formula, 2 means global acromat bend, 1 means double acromat bend, 0 means spreader OFF

flagLH = {'noLH','LH'};
flagIBS = {'noIBS','IBS'};
flagSpreader = {'woS','wS','wS','wS'};
%%
%%%%%%%%%%%%%%%%%%%%%%
% Input Parameters   %
%%%%%%%%%%%%%%%%%%%%%%

%% General
Q = 250e-12;              % total bunch charge in [C]
I0 = 40;                  % initial peak current in [A]
E0 = 6.65e6;              % initial mean energy in [eV]
enx = 0.4e-6;               % normalized horizontal emittance rms in [m rad]
eny = 0.4e-6;               % normalized vertical emittance rms in [m rad]
sigdE = 1;              % natural initial uncorrelated energy spread rms in [eV]
sigd0 = sigdE/E0;  
%%
options = optimset('MaxFunEvals',1000,'TolFun',1e-15,'TolX',1e-15);
% lb = c*Q/(5.5*I0);                                           % RMS bunch length
% 

% D = re^2*(Q/e)/(8*lb*enx);  % constant factor in Bane's approximation