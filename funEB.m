function f = funEB(x,k,b1,tt, params, useDiscretePropDiv, useDiscriteCashDiv)

persistent den eb1 yp psi sqtau sqpi z eta brP et sqt ss ...
    T K cp r q sig beta alpha rhoP1 rP1 exDiscrCashDates cashAmouns 

N = length(tt);
if k == 2 && isempty(sqpi)
    T = params.T;
    K = params.K;
    cp = params.cp;
    r = vertcat(params.r);
    q = vertcat(params.q);
    sig = vertcat(params.sig);
    beta = vertcat(params.beta);
    alpha = vertcat(params.alpha);
    rhoP1 = params.rhoP1;
    rP1 = params.rP1;
    exDiscrCashDates = params.exDiscrCashDates;
    cashAmouns = params.cashAmouns;

    sqpi = sqrt(pi); 
    den(1:N) = 0;
    eb1(1:N) = 1;
    G1(1:N) = 0;
    z(1:N) = 0;
    ss(1:N) = 0;
    eta(1:N) = -rhoP1;
    psi(1:N) = -1;
    sqtau = sqrt(tt);    
    brP(1:N) = rP1;
    rhoP(1:N) = rhoP1;
    yp(1:N) = 0;
    et = exp(tt);
    sqt = sqrt(tt);
end

s = k-1;
tau = tt(k); ss = tt(1:s); dd = tt(k)-ss;
dx = x - b1(1:s); 
gauss = exp(- dx.^2./(4.*dd) );
der = dx./dd;
b1(k) = x;  eb1(k) = exp(x);
psi(k) = - eb1(k).*et(k).*(1 + erf((2*tt(k) + x)./(2*sqt(k))));
z(k) = beta(k).*(alpha(k) - eb1(k));
den = sqrt(pi.*dd);
brP(k) = 2*r(k)/sig(k)^2; rhoP(k) = 2*(r(k)-q(k))/sig(k)^2;

if k < 3    
    yp(2) = (b1(2) - b1(1))./(tt(2) - tt(1));
elseif k < 4
    yp(3) = (3.*b1(3) - 4.*b1(2) + b1(1))./(2.*(tt(3)-tt(2))); 
else
    yp(k) = (11.*b1(k) - 18.*b1(k-1) + 9.*b1(k-2) - 2.*b1(k-3))./(6.*(tt(k)-tt(k-1))); 
end
eta(k) = -(2.*tt(k) + brP(k)).*z(k) - beta(k).*eb1(k).*rhoP(k);

J10 = 0;
if useDiscretePropDiv
    a0 = 2.*ss.*(3 + yp(1:s)) + brP(1:s) - rhoP(1:s);
    a1 = 2.*(-6.*ss.^2 -1 + ss.*(2 + rhoP(1:s) - brP(1:s) + yp(1:s)));
    a2 = -(1 + 4.*ss.^2.*(3 + yp(1:s)));
    a3 = 2.*ss.*(4.*ss.^2 + 1);
    b0 = -2.*ss.*alpha(1:s).*yp(1:s);
    b1 = 2.*alpha(1:s).*(1 + 6.*ss.^2 + ss.*brP(1:s));
    b2 = 4.*ss.^2.*alpha(1:s).*yp(1:s);
    b3 = -2.*ss.*alpha(1:s).*(4.*ss.^2 + 1);
    k2 = ss + 1./(4.*dd); k22 = k2.*k2; k23 = k2.*k22; sqk2 = sqrt(k2);
   
    A  = 1;
    k1 = A + dx./(2.*dd); k12 = k1.*k1;
    prod = (1 + erf(k1./(2.*sqk2))).*exp(k12./(4.*k2));
    prod(isnan(prod)) = 0;
    J2 = sqpi.*prod./(2.*sqk2);
    J2_1 = 1./(2.*k2) + J2.*k1./(2.*k2);
    J2_2 = k1./(4.*k22) + J2.*(2.*k2 + k12)./(4.*k22);
    J2_3 = (4.*k2 + k12)./(8.*k23) + J2.*k1.*(6.*k2 + k12)./(8.*k23);
    J_10_1 = eb1(1:s).*(a0.*J2 + a1.*J2_1 + a2.*J2_2 + a3.*J2_3);

    A  = 0;
    k1 = A + dx./(2.*dd); k12 = k1.*k1;
    prod = exp(k12./(4.*k2)).*(1 + erf(k1./(2.*sqk2)));
    prod(isnan(prod)) = 0;
    J2 = sqpi.*prod./(2.*sqk2);
    J2_1 = 1./(2.*k2) + J2.*k1./(2.*k2);
    J2_2 = k1./(4.*k22) + J2.*(2.*k2 + k12)./(4.*k22);
    J2_3 = (4.*k2 + k12)./(8.*k23) + J2.*k1.*(6.*k2 + k12)./(8.*k23);
    J_10_0 = b0.*J2 + b1.*J2_1 + b2.*J2_2 + b3.*J2_3;
    
    J10 = beta(1:s).*(J_10_0 + J_10_1);  
end

divSum = 0;
if useDiscriteCashDiv
    jm = find(exDiscrCashDates < k);
    if ~isempty(jm)        
        for j = jm
            jj = exDiscrCashDates(j);
            gaussD = gauss(jj)./den(jj);
            k2 = tt(jj) + 1./(4.*dd(jj)); sqk2 = sqrt(k2);
            k22 = k2.*k2;

            a0 = 2*tt(jj)*alpha(jj);
            a1 = a0*(1 - 2*tt(jj)*eb1(jj));
            a2 = -4*a0*tt(jj);
            b0 = 2. - 2*tt(jj)*(1 + eb1(jj));
            b1 = -2*tt(jj)*(3 - 2 *tt(jj)*eb1(jj));
            b2 = 4*a0*tt(jj);
            
            A  = 1;
            k1 = -A + dx(jj)./(2.*dd(jj)); k12 = k1.*k1;
            prod = (1 + erf(k1./(2.*sqk2))).*exp(k12./(4.*k2));
            prod(isnan(prod)) = 0;
            J2 = sqpi.*prod./(2.*sqk2);
            J2_1 = 1./(2.*k2) + J2.*k1./(2.*k2);
            J2_2 = k1./(4.*k22) + J2.*(2.*k2 + k12)./(4.*k22);
            J_10_1 = (a0.*J2 - a1.*J2_1 + a2.*J2_2)/eb1(jj);

            A  = 0;
            k1 = -A + dx(jj)./(2.*dd(jj)); k12 = k1.*k1;
            prod = (1 + erf(k1./(2.*sqk2))).*exp(k12./(4.*k2));
            prod(isnan(prod)) = 0;
            J2 = sqpi.*prod./(2.*sqk2);
            J2_1 = 1./(2.*k2) + J2.*k1./(2.*k2);
            J2_2 = k1./(4.*k22) + J2.*(2.*k2 + k12)./(4.*k22);
            J_10_0 = b0.*J2 + b1.*J2_1 + b2.*J2_2;
    
            Lambda = beta(k).*gaussD.*(J_10_0 + J_10_1);
            divSum = divSum + 2.*alpha(jj)./sig(jj).^2.*cashAmouns(j).*Lambda;
        end    
    end
end

ig1 = [ ( (eta(1:s) + J10).*gauss - eta(k))./den, 0];
add1 = sqtau(k)*2*eta(k)/sqpi;

% integ = trapz(tt(1:k), ig1);
integ = 0.5.*sum((ig1(1:s-1) + ig1(2:s)).*(tt(2:s) - tt(1:s-1)));
f = psi(k) - (integ + add1 + divSum);

if isnan(f)
    ME = MException('function value is NaN');
    throw(ME);
end

end