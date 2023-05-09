function [r] = dLdE(z,m,E)
% Luminosity spectra, in eV/s/eV

a = 0.5;% power index, see 2303.06616
E1 = 100;% Lower cut-off, see 1707.04206
E2 = 5.854289225131452e12;% upper cutoff
Ts = 2E5;% 200 KeV, see 1707.04206
syms x
f=@(x) (x.^-a).*exp(-x);% spectra before normalisation
F=@(x) -gammainc(x,1/2,'upper')*gamma(1/2);
Lbol=Luminosity(z,m);
x1 = E1/Ts;
x2=E2/Ts;

Norm = (Lbol/Ts^(a+1))/(F(x2)-F(x1));
% Make sure input E covers [E1, E2]
if min(E) > 1.2*E1
    error('input E too large')
elseif max(E) < 0.9*E2
    error('input E too small')
end

% Normalised dEtot/dt/dE
r = Norm*(E.^(-a)).*exp(-E/Ts);
r = r.*heaviside(E-E1).*heaviside(E2-E);
end