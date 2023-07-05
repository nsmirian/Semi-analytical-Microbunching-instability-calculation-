
classdef IBS    % S matrix and impedence
    methods(Static)
        function sigddrift_ibs=drift(LD, gamma, enx, betax, sigmaz)
            %-------------------------------% IBS %-----------------------------------%
            %   1- Drift
            run('constants.m')
            D=re^2*(Q/e0)/(8*sigmaz*enx)
            qmax_Dr=sqrt(LD*(Q/e0)*re^2/(2*gamma^(3/2)*enx^(3/2)*sigmaz*sqrt(betax)));
            
            % IBS-induced rms relative energy spread with Bane's approx. of the B-M expression
            sigddrift_ibs = sqrt(2*LD*D*log(qmax_Dr*enx/(re*2*sqrt(2)))/...
                (gamma^2*sqrt(enx*betax./gamma)));
            
        end 
            % 2- Linac
            
        function sigdlinac_ibs=linac(Le,enx,betax,gammai, gammaf,sigmaz)
                run('constants.m')
                D=re^2*(Q/e0)/(8*sigmaz*enx)
            
                % Max angle
                qmax_L =sqrt(Le*(Q/e0)*re^2/(2*enx^(3/2)*sigmaz*sqrt(betax)));
            
                % Coulomb Logarithm
                ln1 =log(qmax_L*enx/(re*2*sqrt(2)));
     
                % IBS-induced rms absolute energy spread in [MeV] with Bane's approximation
                % of B-M expression, cumulated through L1
            
                sigdlinac_ibs =sqrt(D/(((gammaf-gammai)/Le)*sqrt(enx*betax))*((1/sqrt(gammai)...
                -1/sqrt(gammaf))*(4*ln1-6) + (log(gammaf^3)/sqrt(gammaf)-log(gammai^3)/sqrt(gammai))));
            
        end
        function [sigd_BC_w, sigdBC_w_ibs]=DisSec(sigd_BC, beta_in, beta_w,theta,...
                enx, gamma,  alpha_in, alpha_w , DL,Lb, C, sigmaz)
            run('constants.m')
           sigmaz_i=sigmaz;
           sigmaz_w=sigmaz/C
            D=re^2*(Q/e0)/(8*sigmaz_i*enx)
            halfL = DL+Lb;
            %e0 = 1.6e-19;
            %sigd_BC     % uncorrelated energy spread at the entrance of BC
            qmax_BC_i = sqrt(halfL*(Q/e0)*re^2/(2*gamma^(3/2)*enx^(3/2)*sigmaz_i*sqrt(beta_in)));
            qmax_BC_w = sqrt(C*halfL*(Q/e0)*re^2/(2*gamma^(3/2)*enx^(3/2)*sigmaz_w*sqrt(beta_w)));
            
            % optics functions in the first half (_i) and second half (_w) of BC
            etap1 = theta;
            eta1  = theta*(DL+Lb);
            gamma1_in = (1+alpha_in^2)/beta_in;
            gamma1_w  = (1+alpha_w^2)/beta_w;
            
            if switch_bane1 == 1
                
                if theta ~= 0
                    
                    % IBS-induced fractional rms energy spread with Bane's formula for the first half of BC
                    dispH1_i = beta_in*etap1^2+2*alpha_in*eta1*etap1+gamma1_in*eta1^2;
                    aBC1 = vpa(2*D/(gamma^(3/2)*sqrt(enx*beta_in)));
                    bBC1 = vpa(qmax_BC_i*enx/(2*sqrt(2)*re));
                    hBC1 = vpa(dispH1_i*gamma/enx);
                    solBC1 = @(y) vpa(abs( -ei(3*log(sqrt(hBC1*sigd_BC^2 + 1)/bBC1)) + ei(3*log(sqrt(hBC1*y^2 + 1)/bBC1)) + vpa(halfL)*aBC1*hBC1/(2*bBC1^3)));
                    yBC1 = vpa(fminsearch(solBC1,1e-6,options));
                    sigdBC_i_ibs = double(sqrt(yBC1^2 - sigd_BC^2));
                    
                    sigd_BC_w = sqrt(sigd_BC^2 + sigdBC_i_ibs^2);     % Uncorrelated energy spread at the waist of BC
                    
                    % IBS-induced fractional rms energy spread with Bane's formula for the second half of BC
                    dispH1_w = beta_w*etap1^2+2*alpha_w*eta1*etap1+gamma1_w*eta1^2;
                    aBC2 = vpa(2*C*D/(gamma^(3/2)*sqrt(enx*beta_w)));
                    bBC2 = vpa(qmax_BC_w*enx/(2*sqrt(2)*re));
                    hBC2 = vpa(dispH1_w*gamma/enx);
                    solBC2 = @(y) vpa(abs( -ei(3*log(sqrt(hBC2*sigd_BC_w^2 + 1)/bBC2)) + ei(3*log(sqrt(hBC2*y^2 + 1)/bBC2)) + vpa(halfL)*aBC2*hBC2/(2*bBC2^3)));
                    yBC2 = vpa(fminsearch(solBC2,1e-5,options));
                    sigdBC_w_ibs = double(sqrt(yBC2^2 - sigd_BC_w^2));
                    
                elseif theta == 0
                    
                    sigdBC_i_ibs = sqrt(2*halfL*D*log(qmax_BC_i*enx/(re*2*sqrt(2)))./(gamma.^2*sqrt(enx*beta_in./gamma)));
                    sigd_BC_w = sqrt(sigd_BC^2 + sigdBC_i_ibs^2);
                    sigdBC_w_ibs = sqrt(2*halfL*C*D*log(qmax_BC_w*enx/(re*2*sqrt(2)))./(gamma.^2*sqrt(enx*beta_w./gamma)));
                    
                end
                
            elseif switch_bane1 == 0
                
                sigdBC_i_ibs = 0;
                sigdBC_w_ibs = 0;
                sigd_BC_w = sqrt(sigd_BC^2);
                
            end
            
            % IBS-induced energy spread in BC, in [MeV]
            % sigEBC_i = sigdBC_i_ibs*Ef1;
            % sigEBC_w = sigdBC_w_ibs*Ef1;
        end
    end
end
