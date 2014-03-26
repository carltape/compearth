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

            disp([num2str(nn) '/' num2str(niter) ' : prior, current, target:']);
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

            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
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

            % update model
            b = Ga*p;
            mu = g'*icprior*p / (p'*icprior*p + b'*icobs*b);  % Eq. 6.333
            mnew = m - mu*p;
            gold = g;
            pold = p;
            lold = l;

            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
        end

    case 5      % conjugate gradient (polynomial line search)

        for nn = 1:niter
            disp([' iteration ' num2str(nn) ' out of ' num2str(niter) ]);
            m     = mnew;
            Sval  = S(m,dobs,mprior,icobs,icprior);   % misfit value

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
            mu_test   = -2*Sval / sum( g'*icprior*p );
            m_test    = m + mu_test*p;
            Sval_test = S(m_test,dobs,mprior,icobs,icprior);

            % end iteration if the test model is unrealistic
            if ~isreal(Sval_test)
                disp(' polynomial step is TOO FAR');
                S_vec(nn+1:end) = S_vec(nn);
                break;
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

            % get the analytical minimum of the parabola
            if Pa ~= 0, mu = -Pb / (2*Pa); else error('check the input polynomial'); end
            %--------------------

            % update model
            mnew = m + mu*p;
            gold = g;
            pold = p;
            lold = l;

            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
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
            
            disp(sprintf('%i/%i : prior, current, target:',nn,niter));
            disp([mprior mnew mtarget]);
            Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
            Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
            S_vec(nn+1) = S(mnew,dobs,mprior,icobs,icprior);
        end

        % estimated posterior covariance matrix from F_hat, Eq. 6.362
        % compare with cpost = inv(Gpost'*icobs*Gpost + icprior)
        Fhat = F * cprior

end  % case imethod

%==========================================================================
