% Κλείσιμο όλων των ανοιχτών figures και καθαρισμός του workspace, ώστε να
% σιγουρευτεί ότι καμία προηγούμενη εκτέλεση κώδικα δεν θα επηρεάσει την
% επόμενη. Επιπλέον, καθαρισμός του command window για ευκολότερη ανάγνωσή
% οποιουδήποτε μηνύματος τυπωθεί εκεί από το πρόγραμμα
close all;
clear;
clc;

% Ορισμός παραμέτρων σύμφωνα με την εκφώνηση
wvl = 1;
k = 2*pi / wvl;
d = wvl / 2;

% Ορισμός ρυθμίσεων για τον γενετικό αλγόριθμο, σύμφωνα με την εκφώνηση
options = optimoptions('ga', PopulationSize=200, MaxGenerations=1000, ...
    MaxStallGenerations=500, MaxStallTime=200, PlotFcn='gaplotbestf');

[I,~] = ga(@nuAFsum,5,[],[],[],[],[1 1 1 1 1],[4 4 4 4 4],[],[],options);

disp(I);

% Δημιουργία του κατακόρυφου διαγράμματος ακτινοβολίας και εύρεση του
% λόγου του ύψους του 1ου (και υψηλότερου) πλευρικού λοβού σε σχέση με τον
% κύριο, σε καθαρό αριθμό και σε dB
theta = 0:180;
thetarad = deg2rad(theta);
psi = k*d*cos(thetarad);
AF = nuAF(psi, I);
AF = AF/max(AF);
figure;
polarplot(thetarad, AF);
title('Radiation pattern — Vertical plane', FontSize=15, ...
    FontWeight='bold');
set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
figure;
plot(theta, AF);
title('Radiation pattern — Vertical plane', FontSize=15, ...
    FontWeight='bold');
xlabel('θ (°)', FontSize=15, FontWeight='bold');
ylabel('|A|', FontSize=15, FontWeight='bold', Rotation=0);
fslr = max(AF(1:74))/max(AF);
fslrdb = 20*log10(fslr);
disp(fslr);
disp(fslrdb);

% Συνάρτηση υπολογισμού του μέτρου του παράγοντα της μη-ομοιόμορφης
% στοιχειοκεραίας, σύμφωνα με την θεωρητική ανάλυση
function A = nuAF(psi, I)
    a1 = I(1)*(1+exp(1i*9*psi));
    a2 = I(2)*(exp(1i*psi)+exp(1i*8*psi));
    a3 = I(3)*(exp(1i*2*psi)+exp(1i*7*psi));
    a4 = I(4)*(exp(1i*3*psi)+exp(1i*6*psi));
    a5 = I(5)*(exp(1i*4*psi)+exp(1i*5*psi));
    A = abs(a1+a2+a3+a4+a5);
end

