function b1 = computeEB(tau,exDiscrCashDates,cashAmouns)

global T K cp r q sig beta alpha rhoP1 rP1
global useDiscretePropDiv useDiscriteCashDiv

params = struct('T', T, 'K', K, 'cp', cp, 'r', r, 'q', q, 'sig', sig, 'beta', ...
    beta, 'alpha', alpha, 'rhoP1', rhoP1, 'rP1', rP1, ...
    'exDiscrCashDates', exDiscrCashDates, 'cashAmouns', cashAmouns);

N = length(tau);
% options = optimset('fzero');
% options.TolX = 1.e-6;

SbT = K;
b1(1:N) =  log(alpha(1).*SbT./K);
for k = 2:N
    try
        if useDiscretePropDiv || useDiscriteCashDiv
            [b,fval,exitflag,output] = fzero(@(x) funEB(x,k,b1,tau,params, ...
                useDiscretePropDiv,useDiscriteCashDiv), b1(k-1)-0.001);  
        else
            [b,fval, it] = myFzero(k,b1(k-1)-1,b1(k-1)-0.001);  
        end
        b1(k) = b; 
    catch ME
        fprintf('Return by exception %s\n', ME.message)
        return
    end        
end

    function [b, f, it] = myFzero(k, bl, br)
        maxIter = 100;
        tol = 1.e-6;
        fl = funEB(bl,k,b1,tau,params);
        fr = funEB(br,k,b1,tau,params);

        for it=1:maxIter
            if abs(fl) < tol 
                f = fl; b = bl; return
            elseif abs(fr) < tol
                f = fr; b = br; return
            else
                b = 0.5*(bl + br);
                f = funEB(b,k,b1,tau,params);
                if sign(f) == sign(fr)
                    br = b; fr = f;                    
                else
                    bl = b; fl = f;
                end
                if (abs(bl-br) < tol) || (abs(f) < tol) || (it > maxIter)
                    return
                end
            end
        end
    end

end
