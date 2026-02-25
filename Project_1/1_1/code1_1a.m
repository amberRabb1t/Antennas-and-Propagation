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
phi = 0:0.001:2*pi;
theta = 0:0.001:2*pi;

% Ορισμός των τιμών που εξετάζουμε, αριθμητικά για τις πράξεις και σε
% string για την σωστή εμφάνισή τους με σύμβολα (και όχι νούμερα) στα plots
d = [wvl/2 3*wvl/4];
ddisp = ["λ/2" "3λ/4"];

phival = [pi/2 pi/3 pi/6 0];
angledisp = ["0°" "30°" "60°" "90°"];
nphv = length(phival);

dx = [0 -pi/2 -sqrt(3)*pi/2 -pi; 0 -3*pi/4 -3*sqrt(3)*pi/4 -3*pi/2];

% Δημιουργία των οριζόντιων και κατακόρυφων διαγραμμάτων ακτινοβολίας για
% όλους τους συνδυασμούς d και θ, κάνοντας χρήση (σε λούπα που διατρέχει
% όλες τις ζητούμενες τιμές) των κατάλληλων συναρτήσεων που ορίζονται στο
% τέλος του προγράμματος
for n=1:length(d)
    figure;
    for m=1:nphv
        subplot(2, 2, m);
        sgtitle('Radiation pattern — Horizontal plane', FontSize=15, ...
            FontWeight='bold');
        polarplot(Hor(k, d(n), phi, dx(n,m)));
        title(sprintf('d=%s, θ=%s', ddisp(n), angledisp(m)), FontSize=14);
    end
    figure;
    for m=1:nphv
        subplot(2, 2, m);
        sgtitle('Radiation pattern — Vertical plane', FontSize=15, ...
            FontWeight='bold');
        polarplot(Ver(k, d(n), theta, dx(n,m), phival(m)));
        title(sprintf('d=%s, θ=%s', ddisp(n), angledisp(m)), FontSize=14);
        set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    end
end

% Συνάρτηση υπολογισμού του οριζόντιου διαγράμματος ακτινοβολίας, σύμφωνα
% με την θεωρητική ανάλυση
function E = Hor(k, d, phi, dx)
    E = abs(sin(8*k*d*cos(phi)+8*dx)./sin(k*d*cos(phi)/2+dx/2));
end

% Συνάρτηση υπολογισμού του κατακόρυφου διαγράμματος ακτινοβολίας, σύμφωνα
% με την θεωρητική ανάλυση
function E = Ver(k, d, theta, dx, phival)
    E = abs(cos(pi*cos(theta)/2)./sin(theta)).*abs(sin(8*k*d*cos(phival)...
        *sin(theta)+8*dx)./sin(k*d*cos(phival)*sin(theta)/2+dx/2)).* ...
            abs(sin(6*k*d*cos(theta))./sin(k*d*cos(theta)/2));
end

