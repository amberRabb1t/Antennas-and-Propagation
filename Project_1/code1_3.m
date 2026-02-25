% Κλείσιμο όλων των ανοιχτών figures και καθαρισμός του workspace, ώστε να
% σιγουρευτεί ότι καμία προηγούμενη εκτέλεση κώδικα δεν θα επηρεάσει την
% επόμενη. Επιπλέον, καθαρισμός του command window για ευκολότερη ανάγνωσή
% οποιουδήποτε μηνύματος τυπωθεί εκεί από το πρόγραμμα
close all;
clear;
clc;

% Ορισμός απαραίτητων παραμέτρων
c = 2.998*10^8;
f = 10^9;
wvl = c / f;
k = 2*pi / wvl;
l = wvl / 2;

% Υπολογισμός της αυτοαντίστασης κεραίας λ/2 χρησιμοποιώντας τους
% θεωρητικούς τύπους
kl = k*l;
g = double(eulergamma);
a = wvl*10^(-5);
Rm = 60*(g+log(kl)-cosint(kl))+30*sin(kl)*(sinint(2*kl)-2*sinint(kl))+ ...
    30*cos(kl)*(g+log(kl/2)+cosint(2*kl)-2*cosint(kl));
Xm = 60*sinint(kl)+30*cos(kl)*(2*sinint(kl)-sinint(2*kl))-30*sin(kl)* ...
    (2*cosint(kl)-cosint(2*kl)-cosint(2*k*a^2/l));
Z22 = Rm + 1i*Xm;

% α) Δημιουργία του γραφήματος της αμοιβαίας σύνθετης αντίστασης δύο
% παράλληλων διπόλων λ/2 σε απόσταση d χρησιμοποιώντας την συνάρτηση zmpd,
% η οποία ορίζεται στο τελος του προγράμματος
vec = 0:0.001:3*wvl;
Z21m = zmpd(k,vec,l);
vec = vec / wvl;
plot(vec, real(Z21m));
hold on;
plot(vec, imag(Z21m));
yline(0,'--');
title('Z_2_1_m vs d/λ', FontSize=15);
ylabel('Mutual impedance Z_2_1_m (ohms)', FontSize=14, FontWeight='bold');
xlabel('Distance over wavelength d/λ', FontSize=14, FontWeight='bold');
legend('R_2_1_m','X_2_1_m');

% β) Υπολογισμός της σύνθετης αντίστασης εισόδου και σχεδιασμός του
% οριζόντιου διαγράμματος ακτινοβολίας της κεραίας, σύμφωνα με την
% θεωρητική ανάλυση
phi = 0:0.001:2*pi;
d = wvl / 4;
Z21m = zmpd(k,d,l);
Z13m = zmpd(k,2*d,l);
i1overi2 = -Z21m/(Z22+Z13m);
Zin2 = 2*Z21m*i1overi2 + Z22;
fprintf('The input impedance is %g + j%g.\n', real(Zin2), imag(Zin2));
figure;
polarplot(abs(1+2*i1overi2*cos(k*d*cos(phi))));
title(sprintf('Radiation pattern — Horizontal plane\nd=λ/4'), FontSize=14);

% γ) Δημιουργία του γραφήματος του μέτρου του συντελεστή ανάκλασης στην
% είσοδο της κεραίας συναρτήσει της απόστασης d των στοιχείων και
% εντοπισμός της περιοχής τιμών για την οποία είναι μικρότερο του 0.3
vec = 0:0.001:wvl;
Z21m = zmpd(k,vec,l);
Z13m = zmpd(k,2*vec,l);
i1overi2 = -Z21m./(Z22+Z13m);
Zin2 = 2*Z21m.*i1overi2 + Z22;
S11 = abs((Zin2-50)./(Zin2+50));
figure;
plot(vec/wvl, S11);
title('|S_1_1| vs d/λ', FontSize=15);
ylabel('Magnitude of Reflection Coefficient |S_1_1|', FontSize=14, ...
    FontWeight='bold');
xlabel('Distance over wavelength d/λ', FontSize=14, FontWeight='bold');
yline(0.3,'-','0.3',LabelHorizontalAlignment='center',...
    LabelVerticalAlignment='top');
xline(0.457,'-','0.457',LabelOrientation='horizontal', ...
    LabelHorizontalAlignment='center',LabelVerticalAlignment='bottom');
xline(0.607,'-','0.607',LabelOrientation='horizontal', ...
    LabelHorizontalAlignment='center',LabelVerticalAlignment='bottom');

% δ) Υπολογισμός της αντίστασης εισόδου της κεραίας μετά την προσθήκη
% ανακλαστήρα, σύμφωνα με την θεωρητική ανάλυση. Στη συνέχεια, δημιουργία 
% 3D και 2D γραφήματος του μέτρου του συντελεστή ανάκλασης στην είσοδο της
% κεραίας συναρτήσει της απόστασης d των στοιχείων και της απόστασης h του
% ανακλαστήρα από τα αυτά
vec = 0:0.0025:wvl;
n = length(vec);
S11 = zeros(n,n);
Z21m = zmpd(k,vec,l);
Z13m = zmpd(k,2*vec,l);
for h=1:n
    d1 = sqrt((2*vec(h))^2+vec.^2);
    d2 = sqrt((2*vec(h))^2+4*vec.^2);
    Z14m = zmpd(k,2*vec(h),l);
    Z15m = zmpd(k,d1,l);
    Z16m = zmpd(k,d2,l);
    i1overi2 = (Z15m-Z21m)./(Z22+Z13m-Z14m-Z16m);
    Zin2 = 2*Z21m.*i1overi2 + Z22 - 2*Z15m.*i1overi2 - Z14m;
    S11(:,h) = abs((Zin2-50)./(Zin2+50));
end
figure;
surf(vec/wvl,vec/wvl,S11');
title(sprintf(['Magnitude of Reflection Coefficient |S_1_1| for a given'...
    ' d and h']), FontSize=14);
xlabel('d/λ', FontSize=15, FontWeight='bold');
ylabel('h/λ', FontSize=15, FontWeight='bold');
zlabel('|S_1_1|', FontSize=15, FontWeight='bold');
figure;
contour(vec/wvl,vec/wvl,S11',0:0.1:0.3,ShowText='on');
title('Magnitude of Reflection Coefficient |S_1_1| — Contour plot', ...
    FontSize=15);
ylabel('h/λ', FontSize=14, FontWeight='bold');
xlabel('d/λ', FontSize=14, FontWeight='bold');

% Συνάρτηση υπολογισμού της αμοιβαίας σύνθετης αντίστασης δύο
% παράλληλων διπόλων μήκους l σε απόσταση d, χρησιμοποιώντας τους
% θεωρητικούς τύπους
function Z21m = zmpd(k, d, l)
    u0 = k*d;
    u1 = k*(sqrt(d.^2+l^2)+l);
    u2 = k*(sqrt(d.^2+l^2)-l);
    R21m = 120*pi/(4*pi)*(2*cosint(u0)-cosint(u1)-cosint(u2));
    X21m = -120*pi/(4*pi)*(2*sinint(u0)-sinint(u1)-sinint(u2));
    Z21m = R21m + 1i*X21m;
end

