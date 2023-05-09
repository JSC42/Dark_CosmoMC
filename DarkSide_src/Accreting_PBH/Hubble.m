function [r]=Hubble(z)
OmM = 0.3111;
OmR = 9.065340263942517E-5;
OmL = 1 - OmM - OmR;
h = 0.6766;
zp = 1 + z;
H0 = 3.240755744239557e-18 * h;
r = H0 * sqrt(OmL + OmM * (zp.^3) + OmR * (zp.^4));
end