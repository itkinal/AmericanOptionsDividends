%%  American BS model via integral equations
clear all;
close all;
clc;
warning('off', 'all');
set(0, 'DefaultFigureWindowStyle', 'docked');
% profile on
% clearAllMemoizedCaches

global Ar Br Cr Aq Bq Cq As Bs Cs
global T K cp
global r q sig beta alpha rhoP1 rP1
global useDiscretePropDiv useDiscriteCashDiv

loyolagray = 1/255*[169,169,169];
plotInitData = 0;
compEB = 1;

K = 100; S = 100; T = 0.25; cp = '';
N = 50;
tt = linspace(0,T, N);

useDiscretePropDiv = false;
useDiscriteCashDiv = true;

%% Initial values of time-dependent parameters 
Ar = 0.01; Br = 1.; Cr = 0.01;
rr = @(t) Ar*exp(-Br*t) + Cr;  
As = 0.6; Bs = 2; Cs = 0.0;
ssig = @(t) As*exp(-Bs*t) + Cs;  
Aq = 0.02; Bq = 0.5; Cq = -0.01;   
rq = @(t) Aq*exp(-Bq*t) + Cq;  

[exDiscrPropDates, propAmouns, exDiscrCashDates, cashAmouns] = deal([]);
r = flip(rr(tt)); 
sig = flip(ssig(tt));
if useDiscretePropDiv
    exDiscrPropDates = [0.07, 0.12, 0.17, 0.22];
    propAmouns = [0.05,0.04,0.03, 0.02];
end
q = flip(dividends(tt,rq, exDiscrPropDates,propAmouns, exDiscrCashDates,cashAmouns)); 

if useDiscriteCashDiv
    exDiscrCashDates = [0.07, 0.12, 0.17, 0.22];
    cashAmouns = [0.05,0.04,0.03, 0.02].*K;
    for i = 1:length(exDiscrCashDates)
        exDiscrCashDates(i) = 0.5.*integral(@(s) ssig(s).^2, exDiscrCashDates(i), T);
    end
    exDiscrCashDates = flip(exDiscrCashDates);    
    cashAmouns = flip(cashAmouns);
end   

if plotInitData
    plotData(tt,r,q,sig);
end

tau(1:N) = 0; brq(1:N) = 0; br(1:N) = 0; 
for i=1:N
    tau(i) = 0.5.*integral(@(s) ssig(s).^2, tt(i), T);
    bRho(i) = integral(@(s) rr(s) - dividends(s,rq, exDiscrPropDates,propAmouns, [],[]), tt(i), T);
    br(i) = integral(@(s) rr(s), tt(i), T);
end    
tau = flip(tau); tau(1) = 0;
for i=1:length(exDiscrCashDates)
    [c,index] = min(abs(tau-exDiscrCashDates(i)));
    exDiscrCashDates(i) = index;
end
alpha = exp(-tau + flip(bRho));
brr = exp(flip(br));
rhoP1 = r(1)-q(1);
rP1 = r(1);
beta = brr./alpha;

%% Main block, computation of the early exercise boundary
use2int = 0; 
if compEB
    %%%%%%%%%%%%%%%% GIT %%%%%%%%%%%%%%%%%%%%%%
    tic
    b1 = computeEB(tau,exDiscrCashDates,cashAmouns);
    fprintf('\nDone computing EB GIT, elapsed time=%.2g\n', toc);    
    SB = flip(exp(b1).*K./alpha);
    N = length(SB);
    t2 = flip(tt(end:-1:1));

    % if use2int
    %     tic
    %     b1 = computeEB(tau,1);
    %     fprintf('\nDone computing EB GIT, elapsed time=%.2g\n', toc);    
    %     SB1 = flip(exp(b1).*K./alpha);
    % end        

    %%%%%%%%%%%%%%% Binomial tree %%%%%%%%%%%%
    br = integral(@(s) rr(s), tt(1), T)./T; 
    bq = integral(@(s) dividends(s,rq, exDiscrPropDates,propAmouns, [],[]), tt(1), T)./T; 
    bsig = sqrt(integral(@(s) ssig(s).^2, tt(1), T)/T);

    tic
    % last two parameters should be adjusted appropriately to get stable results
    [S1,t,V] = AmericanOption(K,T,br,bq,bsig, 'put',500, 2000);
    Sf = FreeBoundary(S1,t,V,K,'put');
    fprintf('\nDone computing EB BinTree, elapsed time=%.2g\n', toc);    
    FreeBoundaryIndices = [1, find(abs(diff(Sf)) > 1e-5) + 1];   
    Sf1 = interp1(t(FreeBoundaryIndices), Sf(FreeBoundaryIndices), t2, "pchip");

    % MC
    % M = 50000;
    % B = Boundary( K,br,bq,bsig,T,M);
    % tMC = linspace(0,T,M+1);
    % SBMC = interp1(tMC, B', t2, 'spline');

    % tic
    % [optionPrice, SFD] = americanPutFD(100, 100, 0.001, 1, 0.25, 200, 1, 0.01);
    % Sf2 = interp1(linspace(0,T, length(SFD)), SFD, t2, "spline");   
    % fprintf('\nDone computing EB FD, elapsed time=%.2g\n', toc);    

    plotResults(t2, SB, t2, Sf1, t2, (SB-Sf1)./SB);
    save('savedEB.mat', 'b1')
else
    load = matfile('savedEB.mat');
    % do something
end
% profile viewer
