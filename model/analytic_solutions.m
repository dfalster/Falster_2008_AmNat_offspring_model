
%some solution to offspring size model
%D. Falster dec 2007

%General solution: D= 1,Z!=0
clear
syms R Wm Wr Wt s AE b Wa q z
R = (AE/Wm)*s*(Wt/Wr)^-q*(Wa*Wr/Wt/Wm)^(-b*2/(1+(Wm/Wr)^z))
Wt = (1/(s*AE)*Wa^b*Wr^(1-q))^(1/(b-q))
R=subs(R)
dR_dWm =simple(diff(R, Wm))            %fitness derivative    
d2R_dWm2 =simple(diff(dR_dWm,Wm))        %second derivative    

sel_grad = simple(subs(dR_dWm, Wm, Wr))%evaluate for Wm = Wr
d2R = simple(subs(d2R_dWm2, Wm, Wr))

Wr_E = simple(solve(sel_grad, Wr))      %solve =0 to get Wr*
pretty(Wr_E)      %ESS solution --> note scales with Wa as q/q-1

%reproductive output scaling with Wa
syms o1 B1
pretty(simple(subs(Wr_E, AE, o1*Wa^B1)))    %predicted offspring

%------------------------------------------------------------------------
%check solution at maximum --> is it an ESS??
%numerical analysis
S0 = subs(Wr_E, AE, o1*Wa^B1)    %predicted offspring

S1 = subs(sel_grad, Wr, Wr_E)
S1 = subs(S1, AE, o1*Wa^B1)      %selection gradient

S2 = subs(d2R, Wr, Wr_E)
S2 = subs(S2, AE, o1*Wa^B1)      %second derivative 

subs(S0, {Wa, q,z, b, s, B1, o1}, {10^1, 0.2,1, 0.75, 0.8, 1, 0.5*0.895}) 
subs(S1, {Wa, q,z, b, s, B1, o1}, {10^1, 0.2,1, 0.75, 0.8, 1, 0.5*0.895}) 
subs(S2, {Wa, q,z, b, s, B1, o1}, {10^1, 0.2,1, 0.75, 0.8, 1, 0.5*0.895})


%analytic solution to second deriv - can we simplify??
clear
syms R Wm Wr Wt s AE b Wa q z
R = (AE/Wm)*s*(Wt/Wr)^-q*(Wa*Wr/Wt/Wm)^(-b*2/(1+(Wm/Wr)^z))
%Wt = (1/(s*AE)*Wa^b*Wr^(1-q))^(1/(b-q))
dR_dWm =simple(diff(R, Wm))            %fitness derivative    
d2R_dWm2 =simple(diff(dR_dWm,Wm))        %second derivative    
d2R1 = simple(subs(d2R_dWm2, Wm, Wr))
pretty(d2R1)

%sign of second deriv geivn by long term in brackets
F1 = (8-12*b+4*b^2-6*b*z*log(Wa/Wt)-4*b*z+b^2*z^2*log(Wa/Wt)^2+4*b^2*z*log(Wa/Wt))
simple(subs(F1, log(Wa/Wt), 2*(b-1)/b/z))  %substitute for log(wa/wt) based on scaling of Wt with Wa at ESS

%------------------------------------------------ 
%Check CSS show that second derivative  d2R/dr2 - d2R/dm2  > 0 

clear
syms R Wm Wr Wt s AE b Wa q z o3
R = (AE/Wm)*s*(Wt/Wr)^-q*(Wa*Wr/Wt/Wm)^(-b*2/(1+(Wm/Wr)^z))
dR_dWm =simple(diff(R, Wm))            %fitness derivative    
d2R_dWm2 =simple(diff(dR_dWm,Wm))        %second derivative    
d2R_dm2 = simple(subs(d2R_dWm2, Wm, Wr))
pretty(d2R_dm2)

Wt = (1/(s*AE)*Wa^b*Wr^(1-q))^(1/(b-q))
R = (AE/Wm)*s*(Wt/Wr)^-q*(Wa*Wr/Wt/Wm)^(-b*2/(1+(Wm/Wr)^z))
dR_dWr =simple(diff(R, Wr))            %fitness derivative    
d2R_dWr2 =simple(diff(dR_dWr,Wr))        %second derivative    
d2R_dr2 = simple(subs(d2R_dWr2, Wm, Wr))
pretty(d2R_dr2)

d2R_dr2 = 1/4*(4*z*log(Wa*(1/s/AE*Wa^b*Wr^(1-q))^(-1/(b-q)))*b^2+b^2*z^2*log(Wa*(1/s/AE*Wa^b*Wr^(1-q))^(-1/(b-q)))^2+4*b^2-2*z*log(Wa*(1/s/AE*Wa^b*Wr^(1-q))^(-1/(b-q)))*b-4*b*z*log(Wa*(1/s/AE*Wa^b*Wr^(1-q))^(-1/(b-q)))*q-b*z^2*log(Wa*(1/s/AE*Wa^b*Wr^(1-q))^(-1/(b-q)))^2*q-4*b-4*z*b-4*q*b+4*q+2*z*log(Wa*(1/s/AE*Wa^b*Wr^(1-q))^(-1/(b-q)))*q+4*z)*b*AE*s*((1/s/AE*Wa^b*Wr^(1-q))^(1/(b-q))/Wr)^(-q)*(Wa*(1/s/AE*Wa^b*Wr^(1-q))^(-1/(b-q)))^(-b)/Wr^3/(b-q)
%simplify d2R2 by subtsisituing back for Wt
    %Wt = (1/s/AE*Wa^b*Wr^(1-q))^(1/(b-q))
    %Wt^-1 = (1/s/AE*Wa^b*Wr^(1-q))^(-1/(b-q))   
syms Wt    
d2R_dr2 = 1/4*(4*z*log(Wa*Wt^-1)*b^2+b^2*z^2*log(Wa*Wt^-1)^2+4*b^2-2*z*log(Wa*Wt^-1)*b-4*b*z*log(Wa*Wt^-1)*q-b*z^2*log(Wa*Wt^-1)^2*q-4*b-4*z*b-4*q*b+4*q+2*z*log(Wa*Wt^-1)*q+4*z)*b*AE*s*(Wt/Wr)^(-q)*(Wa*Wt^-1)^(-b)/Wr^3/(b-q)
CSS = simple(d2R_dr2-d2R_dm2)
pretty(CSS)

%sign CSS given by term in brackets, let this = S
S= (-2*b+2*q-z*log(Wa/Wt)*b*q+z*b-2*q*b+2*b^2-z*b*q+z*log(Wa/Wt)*b^2)

%At ESS have 
o3 = (exp(2*(b-q)*(1-b)/b/z/(q-1)))^((1-q)/(b-q))
Wt = o3*Wa
log(Wa/Wt)
S = simple(subs(S))
F1 = simple(subs(S, {b, q}, {0.75, 0.2}))
Z = 0:0.01:20
plot(Z, subs(F1, z,Z), 'color','b')


test = log(Wa/Wt) - 2*(1-b)/b/z
subs(test, {b, q,z}, {0.75, 0.2, 1})

S2= (q-b)*(2-2*b-z*b*log(Wa/Wt))-z*b*(q-1)
S2= (q-b)*(2-2*b-z*b*2*(b-1)/b/z)-z*b*(q-1)
simple(subs(S2, {b, q}, {0.75, 0.2}))

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%General solution: D= 1, Z!=1 and variable mutant benefit
% fraction change in mutant benefit given by G, where G =1 is default model
clear
syms R Wm Wr Wtr Wtm s AE b Wa q z G
Wtm = Wtr*(1+G*(Wm/Wr -1))
R = (AE/Wm)*s*(Wtm/Wm)^-q*(Wa/Wtm)^(-b*2/(1+(Wtm/Wtr)^z))

Wtr = (1/(s*AE)*Wa^b*Wr^(1-q))^(1/(b-q))
R=subs(R)
dR_dWm =simple(diff(R, Wm))            %fitness derivative    
d2R_dWm2 =simple(diff(dR_dWm,Wm))        %second derivative    

sel_grad = simple(subs(dR_dWm, Wm, Wr))%evaluate for Wm = Wr
d2R = simple(subs(d2R_dWm2, Wm, Wr))

Wr_E = simple(solve(sel_grad, Wr))      %solve =0 to get Wr*
pretty(Wr_E)      %ESS solution --> note scales with Wa as q/q-1

%reproductive output scaling with Wa
syms o1 B1
W2 = simple(subs(Wr_E, AE, o1*Wa^B1))    %predicted offspring
W1= simple(subs(W2, G, 1));

clear;
b =0.75; q=0.2; o1 = 0.9/20; B1 =1.0; s=1;
W2 =@(Wa,z,G) exp((2*b-2*q-2*b*q+2.*q.^2+4.*q.*G.*b-2.*q.^2.*G+b.*G.*z.*log(Wa).*q+b.*G.*z.*log(1./s./o1./(Wa.^B1))-2.*b.^2.*G)./b./G./z./(q-1))
W1 =@(Wa,z) exp((2.*b-2.*q+2.*b.*q-2.*b.^2+b.*z.*log(Wa).*q+b.*z.*log(1./s./o1./(Wa.^B1)))./b./z./(q-1))

%compare scaling relationship for different values of G
figure;
W = logspace(-4,4, 50);
loglog(W, W1(W,3))
hold on;
loglog(W, W2(W,3, 20), 'Color', 'red', 'LineStyle','--')
loglog(W, W2(W,3, 0.8), 'Color', 'red', 'LineStyle','--')

%plot in relation to z
figure;
Z=logspace(-2,1,50)
semilogy(Z, W2(1,Z,1)./W1(1,Z))
hold on;
plot(Z, W2(1,Z,0.5)./W1(1,Z))
plot(Z, W2(1,Z,0.1)./W1(1,Z))
plot(Z, W2(1,Z,0.8)./W1(1,Z),'Color', 'red')
plot(Z, W2(1,Z,1.1)./W1(1,Z),'Color', 'red')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%General solution for D > 1 (reserve tissues) 
clear
syms R Wm Wr Wt s AE b Wa q z D k
R = (AE/Wm)*s*D^(-k/Wm)*(Wt/D/Wr)^-q*(Wa*Wr/Wt/Wm)^(-b*2/(1+(Wm/Wr)^z))
Wt = (1/(s*AE)*Wa^b*Wr^(1-q)*D^(k/Wr- q))^(1/(b-q))
R=subs(R)
dR_dWm =simple(diff(R, Wm))            %fitness derivative    
sel_grad = simple(subs(dR_dWm, Wm, Wr))%evaluate for Wm = Wr

Wr_E = simple(solve(sel_grad, Wr))      %solve =0 to get Wr*
pretty(Wr_E)      %ESS solution 

    %reproductive output scaling with Wa
    syms o1 B1
    pretty(simple(subs(Wr_E, AE, o1*Wa^B1)))    %predicted offspring

    %try to simplify 
    syms Y
    subs(Wr_E, k*log(D)*(-2*b+2*q+b*z)*exp((-2*b+2*q-b*z*q*(log(Wa)-log(D))-b*z*log(1/s/AE)+2*b^2-2*b*q)/b/z/(-1+q))/b/z/(-1+q), Y)
    pretty(ans)
    Y = k*log(D)*(-2*b+2*q+b*z)*exp((-2*b+2*q-b*z*q*(log(Wa)-log(D))-b*z*log(1/s/AE)+2*b^2-2*b*q)/b/z/(-1+q))/b/z/(-1+q)
    pretty(simple(Y))
    
%plot effect of D on slope
clear
syms s b Wa q z D k x o1 B1  
WR_E1 = exp((2*b-2*q+b*z*log(Wa)*q+b*z*log(1/s/o1/(Wa^B1))-2*b^2+2*b*q)/b/z/(-1+q));
WR_E2 = exp(-(b*z*lambertw(k*log(D)*(-2*b+2*q+b*z)*exp((-2*b+2*q-b*z*q*(log(Wa)-log(D))-b*z*log(1/s/o1/(Wa^B1))+2*b^2-2*b*q)/b/z/(-1+q))/b/z/(-1+q))-q*b*z*lambertw(k*log(D)*(-2*b+2*q+b*z)*exp((-2*b+2*q-b*z*q*(log(Wa)-log(D))-b*z*log(1/s/o1/(Wa^B1))+2*b^2-2*b*q)/b/z/(-1+q))/b/z/(-1+q))-2*b+2*q-b*z*q*(log(Wa)-log(D))-b*z*log(1/s/o1/(Wa^B1))+2*b^2-2*b*q)/b/z/(-1+q))
L =k*log(D)*(-2*b+2*q+b*z)*exp((-2*b+2*q-b*z*q*(log(Wa)-log(D))-b*z*log(1/s/o1/(Wa^B1))+2*b^2-2*b*q)/b/z/(-1+q))/b/z/(-1+q)


F1 = simple(subs(WR_E1, {Wa, q,z, b, s, B1, o1, D},           {x, 0.2,1, 0.75, 0.01, 0.8, 0.5*0.895, 1})) 
F2 = simple(subs(WR_E2, {k, Wa, q,z, b, s, B1, o1, D}, {10^-5, x, 0.2,1, 0.75, 0.01, 0.8, 0.5*0.895, 2}))
f3 = simple(subs(L, {k, Wa, q,z, b, s, B1, o1, D}, {10^-5, x, 0.2,1, 0.75, 0.01, 0.8, 0.5*0.895, 2}))

X = logspace(-8, 5, 100)
loglog(X, subs(F1, x,X), 'color','r')
hold on;
plot(X, subs(F2, x,X), 'color','b')

subs(F2, x,X) / subs(F1, x,X)   % difference converges on D(q/1-q)

ezplot(lambertw(x), [-1/exp(1), 100])

%to understand this difference lets plot some sutff
subplot(4,1,1);  
loglog(X, subs(F1, x,X), 'color','r')
hold on; plot(X, subs(F2, x,X), 'color','b')
ylabel('ESS offspring size(kg)')

subplot(5,1,2);  
loglog(X, subs(F2, x,X) ./ subs(F1, x,X))
hold on;
plot(X, 2^(0.2/0.8), 'color','r')
ylabel('ratio ESS2 / ESS1')

subplot(5,1,3);
loglog(X, subs(f3, x,X))
ylabel('L')

subplot(5,1,4);
loglog(X, lambertw(subs(f3, x,X)))
ylabel('lambertw(L)')

subplot(5,1,5);
loglog(X, exp(lambertw(subs(f3, x,X))))
ylabel('exp(lambertw(L))')


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%General solution with no competitive asymmetry
%D >1 
clear
syms R Wm Wr Wt s AE b Wa q z D k
z=0
R = (AE/Wm)*s*D^(-k/Wm)*(Wt/D/Wr)^-q*(Wa*Wr/Wt/Wm)^(-b*2/(1+(Wm/Wr)^z))
Wt = (1/(s*AE)*Wa^b*Wr^(1-q)*D^(k/Wr- q))^(1/(b-q))
R=subs(R)
dR_dWm =simple(diff(R, Wm))            %fitness derivative    
sel_grad = simple(subs(dR_dWm, Wm, Wr))%evaluate for Wm = Wr
Wr_E = simple(solve(sel_grad, Wr))      %solve =0 to get Wr*
pretty(Wr_E)      %ESS solution --> note scales with Wa as q/q-1

%when D=1
simple(subs(sel_grad, D, 1))
% note -simplifies to give (b-1)/W0r < 0 --> continual selection for
% smaller offspring size

