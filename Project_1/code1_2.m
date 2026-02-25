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

% Ορισμός εύρους γωνιών
theta = 0:0.001:pi;

% Ορισμός των τιμών που εξετάζουμε, αριθμητικά για τις πράξεις και σε
% string για την σωστή εμφάνισή τους με σύμβολα (και όχι νούμερα) στα plots
h = [wvl/4 wvl/2];
hdisp = ["λ/4" "λ/2"];

th0 = [0 pi/6 pi/3 pi/2];
angledisp = ["0°" "30°" "60°" "90°"];

% Δημιουργία του κατακόρυφου διαγράμματος ακτινοβολίας για όλους τους 
% συνδυασμούς h και θ0, κάνοντας χρήση (σε λούπα που διατρέχει όλες τις 
% ζητούμενες τιμές) της κατάλληλης συνάρτησης που ορίζεται στο
% τέλος του προγράμματος
for n=1:length(h)
    figure;
    for m=1:length(th0)
        subplot(2, 2, m);
        sgtitle('Radiation pattern — Vertical plane', FontSize=15, ...
            FontWeight='bold');
        Eth = calcE(k, theta, 50*wvl, th0(m), h(n));
        Eth = Eth/max(Eth);
        polarplot(theta, Eth);
        title(sprintf('h=%s, θ=%s', hdisp(n), angledisp(m)), FontSize=14);
        set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    end
end

% Δημιουργία της γραφικής παράστασης του |Ε| συναρτήσει του θ0 για h=λ/8
% και θ=90° (οριζόντιο επίπεδο) ώστε να φανεί για ποια θ0 μεγιστοποιείται
% το |Ε| και άρα έχουμε μέγιστη ακτινοβολία στο οριζόντιο επίπεδο
degs = 0:0.0001:180;
figure;
plot(degs,calcE(k,pi/2,100*wvl,deg2rad(degs),wvl/8));
title(sprintf('|E| vs θ_0\n(h=λ/8, θ=90°)'), FontSize=14);
xlabel('Dipole rotation angle θ_0 (°)', FontSize=15, FontWeight='bold');
ylabel('Magnitude of Radiation |E|', FontSize=15, FontWeight='bold');

% Συνάρτηση υπολογισμού του μέτρου του Ε για δεδομένο k, θ, r, θ0 και h
% σύμφωνα με την θεωρητική ανάλυση
function E = calcE(k,theta,r,th0,h)
    r1 = sqrt(r^2+h^2-2*r*h*cos(pi/2-theta));
    r2 = sqrt(r^2+h^2-2*r*h*cos(pi/2+theta));
    th1 = acos(r*cos(theta)./r1)-th0;
    th2 = acos(r*cos(theta)./r2)+th0;
    E1 = 1i*cos(pi*cos(th1)/2)./sin(th1).*exp(-1i*k*r1)./r1;
    E2 = -1i*cos(pi*cos(th2)/2)./sin(th2).*exp(-1i*k*r2)./r2;
    E = abs(E1.*cos(th1-theta-th0)+E2.*cos(th2-theta+th0));
end

