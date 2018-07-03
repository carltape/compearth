%
% optimization_method.m
% called by optimization.m
%
% ii: data index
% kk: model parameter index
% nn: model iteration index
% xx: misc index
%
% Carl Tape, 3/2010
%

% testing (run optimization_method.m)
%mnew = minitial;
%Sd_vec = zeros(niter+1,1); Sm_vec = Sd_vec; S_vec = Sd_vec;
%Sd_vec(1) = Sd_0; Sm_vec(1) = Sm_0; S_vec(1) = S_0;

% display initial model
nn = 0;
disp([' iteration ' num2str(nn) ' out of ' num2str(niter) ]);
disp([num2str(nn) '/' num2str(niter) ' : prior, current, target:']);
disp([mprior mnew mtarget]);

switch imethod
    case 1      % Newton (full Hessian)

        for nn = 1:niter
            disp([' iteration ' num2str(nn) ' out of ' num2str(niter) ]);
            m     = mnew;
            dpred = d(m);
            Ga    = G(m);

            % update the model: Tarantola (2005), Eq 6.319
            % (the line-search parameter is assumed to be nu = 1)
            ghat  = Ga'*icobs*(dpred - dobs) + icprior*(m - mprior);  % gradient
            Hhat1 = icprior + Ga'*icobs*Ga;                           % approximate Hessian
            Hhat2 = zeros(nparm,nparm);
            % The ONLY difference between quasi-Newton and Newton is the
            % Hhat2 term. The iith entry of the residual vector is the
            % weight for the corresponding matrix of partial derivatives (G2).
            % Note that the observations are present in Hhat2 but not in Hhat1.
            dwt = icobs*(dpred-dobs);
            for ii=1:ndata
                Hhat2 = dwt(ii) * G2(m,ii);
            end
            Hhat  = Hhat1 + Hhat2;                                   % full Hessian
            disp('Hhat = Hhat1 + Hhat2:');
            for kk=1:nparm
                disp(sprintf('%8.4f %8.4f %8.4f %8.4f   %8.4f %8.4f %8.4f %8.4f + %8.4f %8.4f %8.4f %8.4f',...
                    Hhat(kk,:),Hhat1(kk,:),Hhat2(kk,:)));
            end
            
            if 1==1
                %mnew  = m - inv(Hhat)*ghat;
                dm    = -Hhat\ghat;
                mnew  = m + dm;
            else
                % equivalent formula: see Tarantola and Valette (1982), Eq. 23-35
                mutemp = (Ga*cprior*Ga' + cobs) \ (dobs-dpred + Ga*(m - mprior));
                mnew  = mprior + cprior*Ga'*mutemp;
            end

            % misfit function for new model
            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);

            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
        end 
    
    case 2      % quasi-Newton

        for nn = 1:niter
            disp([' iteration ' num2str(nn) ' out of ' num2str(niter) ]);
            m     = mnew;
            dpred = d(m);
            Ga    = G(m);

            % update the model: Tarantola (2005), Eq 6.319
            % (the line-search parameter is assumed to be nu = 1)
            Hhat  = icprior + Ga'*icobs*Ga;                           % approximate Hessian
            ghat  = Ga'*icobs*(dpred - dobs) + icprior*(m - mprior);  % gradient

            if 1==1
                %mnew  = m - inv(Hhat)*ghat;
                dm    = -Hhat\ghat;
                mnew  = m + dm;
            else
                % equivalent formula: see Tarantola and Valette (1982), Eq. 23-35
                mutemp = (Ga*cprior*Ga' + cobs) \ (dobs-dpred + Ga*(m - mprior));
                mnew  = mprior + cprior*Ga'*mutemp;
            end

            % misfit function for new model
            % note: book-keeping only -- not used within the algorithm above
            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);

            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
        end

    case 3      % steepest descent

        % NOTE: Use Eq. 6.315 as a preconditioner to convert the 
        %       preconditioned steepest descent to the quasi-Newton method.
        F = F0; % identity
        for nn = 1:niter
            disp([' iteration ' num2str(nn) ' out of ' num2str(niter) ]);
            m = mnew;

            % steepest ascent vector (Eq. 6.307 or 6.312)
            dpred = d(m);
            Ga    = G(m);
            g     = cprior*Ga'*icobs*(dpred - dobs) + (m - mprior);

            % search direction (preconditioned by F) (Eq. 6.311)
            p = F*g;

            % update the model
            b    = Ga*p;
            mu   = g'*icprior*p / (p'*icprior*p + b'*icobs*b);  % Eq. 6.314 (Eq. 6.309 if F = I)
            mnew = m - mu*p;    % Eq 6.297

            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
            
            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
        end

    case 4      % conjugate gradient

        for nn = 1:niter
            disp([' iteration ' num2str(nn) ' out of ' num2str(niter) ]);
            m = mnew;

            % steepest ascent vector (Tarantola, 2005, Eq. 6.312)
            dpred = d(m);
            Ga    = G(m);
            g     = cprior*Ga'*icobs*(dpred - dobs) + (m - mprior);

            % search direction (Tarantola, 2005, Eq. 6.329)
            l = F0*g;
            if nn == 1 
               alpha = 0; p = g;
            else
               alpha = (g-gold)'*icprior*l / (gold'*icprior*lold);  % Eq. 6.331-2
               p = l + alpha*pold;
            end

            % calculate step length
            b = Ga*p;
            mu = g'*icprior*p / (p'*icprior*p + b'*icobs*b);  % Eq. 6.333
            
            % update model
            mnew = m - mu*p;
            gold = g;
            pold = p;
            lold = l;

            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
            
            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
        end

    case 5      % conjugate gradient (polynomial line search)

        for nn = 1:niter
            disp([' iteration ' num2str(nn) ' out of ' num2str(niter) ]);
            m     = mnew;
            Sval  = S(m,dobs,mprior,icobs,icprior);   % misfit value
            
%             if nn > 1
%                 if Sval > S_vec(nn-1)
%                     Sval, S_vec(nn-1)
%                     disp('new model has larger misfit, so stop at previous model');
%                     Sd_vec(nn:end) = Sd_vec(nn-1);
%                     Sm_vec(nn:end) = Sm_vec(nn-1);
%                     S_vec(nn:end)  = S_vec(nn-1);
%                     break
%                 end
%             end
            
            % steepest ascent vector (Tarantola, 2005, Eq. 6.312)
            dpred = d(m);
            Ga    = G(m);
            g     = cprior*Ga'*icobs*(dpred - dobs) + (m - mprior);

            % search direction (Tarantola, 2005, Eq. 6.329)
            l = F0*g;
            if nn == 1 
               alpha = 0; p = -g;
            else
               alpha = (g - gold)'*icprior*l / (gold'*icprior*lold);  % Eq. 6.331-2
               p = l + alpha*pold;
            end
            
            %--------------------
            % test model using quadratic extrapolation: Tape-Liu-Tromp (2007)
            % (compared with the previous CG algorithm, we do not use b = Ga*b to get mu)
            mu_test   = -2*Sval / sum( g'*icprior*p );
            m_test    = m + mu_test*p;
            Sval_test = S(m_test,dobs,mprior,icobs,icprior);

            % end iteration if the test model is unrealistic
            if ~isreal(Sval_test)
                disp('polynomial step is TOO FAR');
                S_vec(nn+1:end) = S_vec(nn);
                break
            end

            % determine coefficients of quadratic polynomial (ax^2 + bx + c),
            % using the two points and one slope
            x1 = 0;
            x2 = mu_test;
            y1 = Sval;
            y2 = Sval_test;
            g1 = sum( g'*icprior*p );
            Pa = ((y2 - y1) - g1*(x2 - x1)) / (x2^2 - x1^2);
            Pb = g1;
            Pc = y1 - Pa*x1^2 - Pb*x1;

            % get the mu value associated with analytical minimum of the parabola
            if Pa ~= 0, mu = -Pb / (2*Pa); else error('check the input polynomial'); end
            %--------------------

            % update model
            mnew = m + mu*p;
            gold = g;
            pold = p;
            lold = l;

            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
            
            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
        end   

    case 6      % variable metric (matrix version)

        F = F0;
        for nn = 1:niter
            disp([' iteration ' num2str(nn) ' out of ' num2str(niter) ]);
            m = mnew;

            % steepest ascent vector
            dpred = d(m);
            Ga    = G(m);
            g     = cprior*Ga'*icobs*(dpred - dobs) + (m - mprior);

            % update the preconditioner F
            if nn > 1
                dg = g - gold;                            % Eq. 6.341
                v  = F*dg;                                % Eq. 6.355
                u  = dm - v;                              % Eq. 6.341
                F  = F + u*u'*icprior / (u'*icprior*dg);  % Eq. 6.356
            end

            % preconditioning search direction (Eq. 6.355)
            p = F*g;

            % update the model
            b    = Ga*p;                                        % Eq. 6.333
            mu   = g'*icprior*p / (p'*icprior*p + b'*icobs*b);  % Eq. 6.333
            dm   = -mu*p;                                       % Eq. 6.355
            mnew = m + dm;
            gold = g; 
            
            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
            
            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
        end

        % estimated posterior covariance matrix from F_hat, Eq. 6.362
        % compare with cpost = inv(Gpost'*icobs*Gpost + icprior)
        Fhat = F * cprior

    case 7      % variable metric (vector version)

        u = NaN(nparm,niter-1);
        v = NaN(nparm,niter-1);
        for nn = 1:niter
            m = mnew;

            % steepest ascent vector
            delta = d(m);
            Ga    = G(m);
            g     = cprior*Ga'*icobs*(delta - dobs) + (m - mprior);

            % update the preconditioner F (Section 6.22.8)
            if nn > 1
                dg = g - gold;
                F_dg = F0*dg;
                % compute u(k) and v(k)
                for jj = 1:(nn-2)  % loop is entered only if nn >= 3
                    vtmp = dg'*icprior*u(:,jj);
                    F_dg = F_dg + vtmp/v(jj) * u(:,jj);     % Eq. 6.347
                end
                u(:,nn-1) = dm - F_dg;                      % Eq. 6.341
                v(nn-1) = dg'*icprior*u(:,nn-1);            % Eq. 6.348
            end

            % preconditioning search direction p = F_g (Eq. 6.340, Eq. 6.347)
            p = F0*g;
            for jj = 1:(nn-1)   % loop is entered only if nn >= 2
                p = p + g'*icprior*u(:,jj)/v(jj) * u(:,jj);
            end

            % update the model
            b    = Ga*p;
            mu   = g'*icprior*p / (p'*icprior*p + b'*icobs*b);  % Eq. 6.333
            dm   = -mu*p;
            mnew = m + dm;
            gold = g; 

            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
        end

        % estimated posterior covariance matrix from F_hat, Eq. 6.362
        F = vm_F(F0,icprior,u,v);
        Fhat = F * cprior

    case 8      % square-root variable metric (matrix version)

        % We should have F = Shat * Shat' * icprior
        F = F0;
        F0hat = F0*cprior;
        S0hat = sqrtm(F0hat);  % works if Fhat is symmetric positive definite
        % check: a0 = rand(2,2); a = a0*a0', b = sqrtm(a), b*b
        Shat = S0hat;
        Fhat = F0hat;

        % store these for checking only
        a_vec  = zeros(niter-1,1);
        b_vec  = zeros(niter-1,1);
        nu_vec = zeros(niter-1,1);
        w_mat  = zeros(nparm,niter-1);

        for nn = 1:niter
            m = mnew;

            % gradient
            delta = d(m);
            Ga    = G(m);
            ghat  = Ga'*icobs*(delta - dobs) + icprior*(m - mprior);

            % update the preconditioner F (using S); see Hull and Tapley (1977)
            if nn >= 2
                dghat = ghat - ghat_old;
                yhat  = mu*ghat_old + dghat;
                w     = Shat'*yhat;
                a     = yhat'*Fhat*dghat;
                b     = w'*w;
                nu    = srvm_nu(a,b);
                Shat  = Shat*(eye(nparm) - nu/a*w*w' );
                Fhat  = Shat*Shat';

                % for checking only
                a_vec(nn-1) = a;
                b_vec(nn-1) = b;
                nu_vec(nn-1) = nu;
                w_mat(:,nn-1) = w;
            end

            % preconditioned gradient
            p = Fhat*ghat;

            % update the model
            c    = Ga*p;
            mu   = ghat'*p / (p'*icprior*p + c'*icobs*c);  % Eq. 6.333
            %mu   = g'*icprior*p / (p'*icprior*p + c'*icobs*c);  % Eq. 6.333
            dm   = -mu*p;
            mnew = m + dm;
            ghat_old = ghat; 

            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
        end

        % estimated posterior covariance matrix from F_hat, Eq. 6.362
        Fhat
        
        % estimated posterior covariance matrix computed from stored vectors and scalars
        Fhat_check = srvm_Fhat(S0hat,niter-1,nu_vec,a_vec,w_mat)
        a = a_vec;
        b = b_vec;
        nu = nu_vec;
        w = w_mat;

    case 9      % square-root variable metric (vector version)

        % We should have F = Shat * Shat' * icprior
        F = F0;
        F0hat = F0*cprior;
        S0hat = sqrtm(F0hat);
        Shat = S0hat;

        % initialize vectors and scalars
        a  = zeros(niter-1,1);
        b  = zeros(niter-1,1);
        nu = zeros(niter-1,1);
        w  = zeros(nparm, niter-1);

        for nn = 1:niter
            m = mnew;

            % gradient
            delta = d(m);
            Ga    = G(m);
            ghat  = Ga'*icobs*(delta - dobs) + icprior*(m - mprior);

            % update the preconditioner F (using S)
            % see Hull and Tapley (1977)
            if nn >= 2
                dghat = ghat - ghat_old;
                yhat  = mu*ghat_old + dghat;

                w(:,nn-1) = srvm_Shat_chi(yhat,nn-2,S0hat,nu,a,w,1);
                beta      = w(:,nn-1) - mu*ShatT_ghat;
                a(nn-1)   = transpose(w(:,nn-1))*beta(:);
                b(nn-1)   = transpose(w(:,nn-1))*w(:,nn-1);
                nu(nn-1)  = srvm_nu(a(nn-1),b(nn-1));
            end

            % update the search direction (does nothing for nn=1)
            ShatT_ghat = srvm_Shat_chi(ghat,nn-1,S0hat,nu,a,w,1);
            p          = srvm_Shat_chi(ShatT_ghat,nn-1,S0hat,nu,a,w,0);

            % update the model
            c        = Ga*p;
            mu       = ghat'*p / (p'*icprior*p + c'*icobs*c);  % Eq. 6.333
            %mu       = g'*icprior*p / (p'*icprior*p + c'*icobs*c);  % Eq. 6.333
            dm       = -mu*p;
            mnew     = m + dm;
            ghat_old = ghat; 

            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
        end

        % estimated posterior covariance matrix computed from stored vectors and scalars
        Fhat = srvm_Fhat(S0hat,niter-1,nu,a,w)          
        
end  % case imethod

%==========================================================================
