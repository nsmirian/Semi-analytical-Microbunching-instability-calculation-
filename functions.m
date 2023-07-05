%functins:

k = @(lambda) 2*pi./lambda;                                  % uncompressed wave number in [1/m]
bk0 = @(lambda) sqrt(2*e0*c./(I0*lambda));                    % initial bunching factor for density modulation, from shot noise
bm0 = @(lambda) 0;                                           % initial bunching factor for energy modulation, from shot noise
bf0 = @(lambda) [bk0(lambda) 0; bm0(lambda) 0];
%options = optimset('MaxFunEvals',1000,'TolFun',1e-15,'TolX',1e-15);


%%
%--------------------------------% LH %-----------------------------------%
%-------------------------------------------------------------------------%
DeltaE_lh=10e3;
sigdE_lh=10e3;
Blh = 1;            % ratio of laser over electron beam radius (>=1) in the LH
% Calculate the hypergeometric function in its integral form. The integral
% runs over the ratio of laser/e-beam radius: assume here from 1 to 10) for
% the definite integral.
A1 = @(lambda, C, R56, Ef) abs((2*pi*C./lambda)*R56*DeltaE_lh/Ef);
J01_LH = @(r,lambda, C, R56, Ef) besselj(0,A1(lambda, C, R56, Ef).*exp(-r.^2/(4*Blh^2)));
J11_LH = @(r,lambda, C, R56, Ef) besselj(1,A1(lambda, C, R56, Ef).*exp(-r.^2/(4*Blh^2)));
S01_LH = @(lambda, C, R56, Ef) integral(@(r)r.*exp(-r.^2/2).*J01_LH(r,lambda, C, R56, Ef),0,Inf);
S11_LH = @(lambda, C, R56, Ef) integral(@(r)r.*exp(-r.^2/2).*exp(-r.^2/(4*Blh^2)).*J11_LH(r,lambda,C, R56, Ef),0,Inf);

if sigdE_lh ~= 0 && switch_lh == 1
    S01_LH = @(lambda, C, R56, Ef) S01_LH(lambda, C, R56, Ef);
    S11_LH = @(lambda, C, R56, Ef) S11_LH(lambda, C, R56, Ef);
elseif sigdE_lh == 0 || switch_lh == 0
    S01_LH = @(lambda, C, R56, Ef) 1;
    S11_LH = @(lambda, C, R56, Ef) 1;
end