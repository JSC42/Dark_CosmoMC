function [r] = dLdE(z,m,E)
% Luminosity spectra, in eV/s/eV

a = 0.5;% power index, see 2303.06616
E1 = 100;% Lower cut-off, see 1707.04206
E2 = 5.854289225131452e12;% upper cutoff
Ts = 2E5;% 200 KeV, see 1707.04206

% Only allow E to be an array
if length(E)<3
    error('E must be an array\n')
end
% Make sure input E covers [E1, E2]
if min(E) > 1.2*E1
    error('input E too large')
elseif max(E) < 0.9*E2
    error('input E too small')
end

f = (E.^(-a)).*exp(-E/Ts);
f = f.*heaviside(E-E1).*heaviside(E2-E);

% Normalise
Lb = Luminosity(z,m);
Norm = Lb/trapz(E,f);
r = Norm*f;

end