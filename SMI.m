


classdef SMI    % S matrix and impedence
    methods(Static)
        %-------------------------------% LSC %-----------------------------------%
        %- 1- DRIFT
        function [S_Drift, Z_LSC_i]=LSCDrift(LD,betax, betay, enx, eny, energy, I)
            %function [S_D, Z_L]=LSCDrift(cell_args)
            run('constants.m')
            run('functions.m')
            gamma=energy/me;
            % Average beam radius in [m] and related quantities
            rb0 = 0.8735*(sqrt(enx*betax/gamma)+sqrt(eny*betay/gamma));
            a0 =  @(lambda) k(lambda)*rb0/gamma;
            I01 = @(lambda) besseli(1,a0(lambda));
            K00 = @(lambda) besselk(0,a0(lambda));
            K01 = @(lambda) besselk(1,a0(lambda));
            aw0 =  @(lambda) k(lambda)*rw0/gamma;
            I00_w = @(lambda) besseli(0,aw0(lambda));
            K00_w = @(lambda) besselk(0,aw0(lambda));
            
            
            % 3-D LSC impedance averaged over transverse dimension (_rr stays for
            % chamber shielding)
            Z_LSC = @(lambda) 1i*(Z0./(pi*k(lambda)*rb0^2)).*(1-2*I01(lambda).*K01(lambda));
            Z_LSC_rr = @(lambda) 1i*(Z0./(pi*k(lambda)*rb0^2)).*(1-2*(I01(lambda)...
                /I00_w(lambda)).*(I00_w(lambda).*K01(lambda)+I01(lambda).*K00_w(lambda)));
            
            % LSC-induced energy modulation amplitude in me unit
            Z_LSC_i = @(lambda) Z_LSC(lambda)*LD;
            %Dgamma = @(lambda, I, LD ) -(4*pi/Z0)*(I/IA)*bk0(lambda).*Z_LSC(lambda)*LD;
            % S-matrix^2
            S_Drift = @(lambda) [1 0;-Z_LSC_i(lambda)*I/energy 1];
            
        end
        %%
        %-----------------------------------------------------------------------
        %%-----------------------------------------------------------------------
        function [S_Li , Z1_LSC_i]=LSCDlinac(L,betax, betay, enx, eny, Ei, Ef, I)
            run('constants.m')
            run('functions.m')
            G1 = (Ef-Ei)/L;                       % L1 average accelerating gradient in [MeV/m]
            gamma1 = @(s1) (Ei+G1*s1)/me;           % Lorentz factor for mean energy along L1
            gammaf = Ef/me;                        % Lorentz factor for mean energy at the END of L1
            gammam = (gammaf - Ei/me)/2;          % Mean Lorentz factor of L1
            
            % Average beam radius in [m] and related quantities
            rb1 = @(s1) 0.8735*(sqrt(enx*betax./gamma1(s1))+sqrt(eny*betay./gamma1(s1)));
            a1 =  @(s1,lambda) k(lambda).*rb1(s1)/gamma1(s1);
            I11 = @(s1,lambda) besseli(1,a1(s1,lambda));
            K10 = @(s1,lambda) besselk(0,a1(s1,lambda));
            K11 = @(s1,lambda) besselk(1,a1(s1,lambda));
            aw1 =  @(s1,lambda) k(lambda)*rw1./gamma1(s1);
            I10_w = @(s1,lambda) besseli(0,aw1(s1,lambda));
            K10_w = @(s1,lambda) besselk(0,aw1(s1,lambda));
            
            % 3-D LSC impedance averaged over transverse dimension (_rr stays for
            % chamber shielding)
            Z1_LSC = @(s1,lambda)...
                1i.*(Z0./(pi*k(lambda).*rb1(s1).^2)).*(1-2*I11(s1,lambda).* K11(s1,lambda));
            Z1_LSC_rr = @(s1,lambda) 1i.*(Z0/(pi*k(lambda).*rb1(s1).^2)).*(1-2*(I11(s1,lambda)./I10_w(s1,lambda)).*(I10_w(s1,lambda).*K11(s1,lambda)+I11(s1,lambda).*K10_w(s1,lambda)));
            
            % LSC-induced energy modulation amplitude in me unit
            %Dgamma1 = @(lambda) -(4*pi/Z0)*(I0/IA).*b0(lambda).*integral(@(s1)Z1_LSC(s1,lambda),0,L1);
            Z1_LSC_real = @(s1,lambda) real(Z1_LSC(s1,lambda));
            Z1_LSC_imag = @(s1,lambda) imag(Z1_LSC(s1,lambda));
            
            % Z_LSC impedance integrated through L0
            Z1_LSC_i = @(lambda) integral(@(s1) Z1_LSC(s1,lambda),0,L);
            Z1_LSC_abs = @(lambda) abs(integral(@(s1) Z1_LSC(s1,lambda),0,L));
            % S-matrix of LinaC
            S_Li = @(lambda) [1 0; -Z1_LSC_i(lambda)*I./Ef Ei/Ef];
            
            
            
        end
        %%  
        function S_BC =LSC_CSRDS(DL,Lb,beta_in,beta_w, enx, eny, Ef, I, theta, C,...
                D12, D13, D14, sigd_BC_w, k1);
            run('constants.m')
            run('functions.m')
            
            halfL = DL+Lb;                                            % half length of the compressor,
            eta_max = theta*(Lb+DL);                                 % maximum dispersion function in the chicane, in [m]
            R56 = -2*theta^2*((2/3)*Lb+DL);                          % R56 of BC01-chicane in [m]
            h = abs((1-1/C)/R56);                                 % linear energy chirp at BC0
            R = Lb/theta;                                             % BC0 dipoles bending radius in [m]                                             % Delay Line chicane bending radius in [m]
            %sigd_BC0_cor = abs(h*Q*c/(I0*sqrt(2*pi)));                  % correlated fractional energy RMS spread at BC0
            
            H01 = (theta^2*beta_in);                                   % horizontal H-function at the exit of first dipole of BC0-DS1, in [m]
            H02 = (theta^2*beta_w+2*eta_max^2/beta_w);               % horizontal H-function at the exit of second+third dipole of BC0-DS1, in [m]
            I1 = I*C; 
            k2 = @(lambda) k1(lambda)*C;
            gamma=Ef/me;
            % Average beam radius in [m] and related quantities
            rbBC_i = 0.8735*(sqrt(enx*beta_in/gamma) + sqrt(eny*beta_in/gamma));
            sBC_i = sqrt(enx*beta_in/gamma) + sqrt(eny*beta_in/gamma);
            spBC_i = sqrt(enx/(gamma*beta_in)) + sqrt(eny/(gamma*beta_in));
            rbBC_w = 0.8735*(sqrt(enx*beta_w/gamma) + sqrt(eny*beta_w/gamma));
            sBC_w = sqrt(enx*beta_w/gamma) + sqrt(eny*beta_w/gamma);
            spBC_w = sqrt(enx/(gamma*beta_w)) + sqrt(eny/(gamma*beta_w));
            aBC_i =  @(lambda) k(lambda)*rbBC_i/gamma;
            aBC_w =  @(lambda) k(lambda)*rbBC_w/gamma;
            I1_BC_i = @(lambda) besseli(1,aBC_i(lambda));
            I1_BC_w = @(lambda) besseli(1,aBC_w(lambda));
            K0_BC_i = @(lambda) besselk(0,aBC_i(lambda));
            K0_BC_w = @(lambda) besselk(0,aBC_w(lambda));
            K1_BC_i = @(lambda) besselk(1,aBC_i(lambda));
            K1_BC_w = @(lambda) besselk(1,aBC_w(lambda));
            awBC =  @(lambda) k(lambda)*rw0/gamma;
            I0_w_BC = @(lambda) besseli(0,awBC(lambda));
            K0_w_BC = @(lambda) besselk(0,awBC(lambda));
            
            % 3-D LSC impedance averaged over transverse dimension (_rr stays for
            % chamber shielding)
            Z_LSC_BC_i = @(lambda) 1i*(Z0./(pi*k(lambda)*rbBC_i^2)).*(1-2*I1_BC_i(lambda).*K1_BC_i(lambda));
            Z_LSC_BC_w = @(lambda) 1i*(Z0./(pi*k(lambda)*rbBC_w^2)).*(1-2*I1_BC_w(lambda).*K1_BC_w(lambda));
            Z_LSC_rr_BC_i = @(lambda) 1i*(Z0./(pi*k(lambda)*rbBC_i^2))*(1-2*(I1_BC_i(lambda)/I0_w_BC(lambda))*(I0_w_BC(lambda)*K1_BC_i(lambda)+I1_BC_i(lambda)*K0_w_BC(lambda)));
            Z_LSC_rr_BC_w = @(lambda) 1i*(Z0./(pi*k(lambda)*rbBC_w^2))*(1-2*(I1_BC_w(lambda)/I0_w_BC(lambda))*(I0_w_BC(lambda)*K1_BC_w(lambda)+I1_BC_w(lambda)*K0_w_BC(lambda)));
            
            
            
            if C ~= 1
                
                % CSR impedance of BC first half
                
                % CSR in dipole-1
                ZCSR_BC_1 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb*k(lambda).^(1/3)/R^(2/3)).*exp(-(k(lambda)*theta*sBC_i).^2);
                % CSR in dipole-2
                ZCSR_BC_2_i = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(C/(1+C))*(Lb*(((2*C)/(1+C))*k1(lambda)).^(1/3)/R^(2/3)).*exp(-(((2*C)/(1+C))*k1(lambda)*DL*theta*spBC_i).^2);
                % CSR from dipole-2
                ZCSR_BC_2_e = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(1/(1+C))*(Lb*((2/(1+C))*(k1(lambda))).^(1/3)/R^(2/3)).*exp(-((2/(1+C))*k1(lambda)*D12*theta*spBC_i).^2);
                
                % CER impedance of BC first half
                
                % CER from dipole-1
                ZCER_BC_1 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(DL,lambda*gamma^2/(2*pi))/(R^(2/3)*lambda.^(1/3)))).*exp(-(k1(lambda)*theta*sBC_i).^2);
                % CER from dipole-2
                ZCER_BC_2_i = @(lambda) switch_cer*(Z0/(2*pi))*(C/(1+C))*log((min(D12,((1+C)/(2*C))*(lambda*gamma^2)/(2*pi))/(R^(2/3)*(lambda*(1+C)/(2*C)).^(1/3)))).*exp(-(((2*C)/(1+C))*k1(lambda)*DL*theta*spBC_i).^2);
                % CER from exit of dipole-2
                ZCER_BC_2_e = @(lambda) switch_cer*(Z0/(2*pi))*(1/(1+C))*log((min(D12,((1+C)/2).*(lambda*gamma^2)/(2*pi))/((R^(2/3).*(lambda*(1+C)/2).^(1/3))))).*exp(-((2/(1+C))*k1(lambda)*D12*theta*spBC_w).^2);
                
                % Total impedance of BC first half
                ZBC_fh = @(lambda) ZCSR_BC_1(lambda) + ZCSR_BC_2_i(lambda) + ZCSR_BC_2_e(lambda) + ZCER_BC_1(lambda) + ZCER_BC_2_i(lambda) + ZCER_BC_2_e(lambda);
                S_BC_fh = @(lambda) [1 0; -ZBC_fh(lambda)*I/Ef 1];
                
                % Amplification & damping functions
                F1_ld = @(lambda) exp(-0.5*(C*R56*k1(lambda)*sigd_BC_w).^2)*S01_LH(lambda, C, R56, Ef);
                G1_ld = @(lambda) exp(-0.5*(C*R56*k1(lambda)*sigd_BC_w).^2)*S11_LH(lambda, C, R56, Ef)*Ef*(C*R56*k1(lambda)*sigd_BC_w^2);
                
                % Compression matrix
                S_BC_c = @(lambda) [F1_ld(lambda) 1i*F1_ld(lambda)*C*R56.*k1(lambda);...
                    1i*G1_ld(lambda)*C/Ef (C*F1_ld(lambda)-C^2*G1_ld(lambda)*R56.*k1(lambda)/Ef)];
                
                % CSR impedance of BC second half
                
                % CSR in dipole-3
                ZCSR_BC_3 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb*k2(lambda).^(1/3)/R^(2/3)).*exp(-(k2(lambda)*theta*sBC_w).^2);
                % CSR in dipole-4
                ZCSR_BC_4 = @(lambda) switch_csr*(Z0/(4*pi))*A_csr*(Lb*k2(lambda).^(1/3)/R^(2/3));
                
                % CER impedance of BC second half
                
                % CER from dipole-3
                ZCER_BC_3 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(D13,lambda*gamma^2/(2*pi))/(R^(2/3)*lambda.^(1/3)))).*exp(-(k2(lambda)*theta*sBC_w).^2);
                % CER from dipole-4
                ZCER_BC_4 = @(lambda) switch_cer*(Z0/(2*pi))*log((min(D14,lambda*gamma^2/(2*pi))/(R^(2/3)*lambda.^(1/3))));
                
                % Total impedance of BC second half
                ZBC_sh = @(lambda) ZCSR_BC_3(lambda) + ZCSR_BC_4(lambda) + ZCER_BC_3(lambda) + ZCER_BC_4(lambda);
                S_BC_sh = @(lambda) [1 0; -ZBC_sh(lambda)*I1/Ef 1];
                
            else
                
                % LSC impedance of BC first half
                ZBC_fh = @(lambda) Z_LSC_BC_i(lambda)*halfL;
                S_BC_fh = @(lambda) [1 0; -ZBC_fh(lambda)*I/Ef 1];
                
                S_BC_c = @(lambda) [1 0; 0 1];
                
                % LSC impedance of BC second half
                ZBC_sh = @(lambda) Z_LSC_BC_w(lambda)*halfL;
                S_BC_sh = @(lambda) [1 0; -ZBC_sh(lambda)*I1/Ef 1];
            end
            % S-matrix of BC0
            S_BC = @(lambda) S_BC_sh(lambda)*S_BC_c(lambda)*S_BC_fh(lambda);
        end
    end
end