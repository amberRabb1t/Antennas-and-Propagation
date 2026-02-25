% Συνάρτηση υπολογισμού του αθροίσματος των τιμών του μέτρου του παράγοντα
% της μη-ομοιόμορφης στοιχειοκεραίας για γωνίες 0-70°, ανά 1°
function A = nuAFsum(p)
    wvl = 1;
    k = 2*pi / wvl;
    d = wvl / 2;
    theta = 0:70;
    theta = deg2rad(theta);
    psi = k*d*cos(theta);
    a1 = p(1)*(1+exp(1i*9*psi));
    a2 = p(2)*(exp(1i*psi)+exp(1i*8*psi));
    a3 = p(3)*(exp(1i*2*psi)+exp(1i*7*psi));
    a4 = p(4)*(exp(1i*3*psi)+exp(1i*6*psi));
    a5 = p(5)*(exp(1i*4*psi)+exp(1i*5*psi));
    A = sum(abs(a1+a2+a3+a4+a5));
end

