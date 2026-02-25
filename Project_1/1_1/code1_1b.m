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

% Ορισμός εύρους γωνιών σε μορφή συντεταγμένων
[phi, theta] = meshgrid(0:0.009:2*pi, 0:0.009:pi);

% Ορισμός των τιμών που εξετάζουμε, αριθμητικά για τις πράξεις και σε
% string για την σωστή εμφάνισή τους με σύμβολα (και όχι νούμερα) στα plots
d = [wvl/2 3*wvl/4];
ddisp = ["λ/2" "3λ/4"];

angledisp = ["0°" "30°" "60°" "90°"];

dx = [0 -pi/2 -sqrt(3)*pi/2 -pi; 0 -3*pi/4 -3*sqrt(3)*pi/4 -3*pi/2];

% Δημιουργία των 3D στερεών ακτινοβολίας για όλους τους συνδυασμούς d και
% θ, κάνοντας χρήση (σε λούπα που διατρέχει όλες τις ζητούμενες τιμές) των
% κατάλληλων συναρτήσεων που ορίζονται στο τέλος του προγράμματος. Για να
% χρησιμοποιηθεί η συνάρτηση surf (η οποία ζωγραφίζει το στερεό), γίνεται
% μετάβαση από σφαιρικές σε καρτεσιανές συντεταγμένες
for axeq=1:2
    for n=1:length(d)
        figure;
        for m=1:length(angledisp)
            subplot(2, 2, m);
            if axeq==2
                sgtitle(sprintf('Radiation pattern — axis equal'), ...
                    FontSize=15, FontWeight='bold');
            else
                sgtitle(sprintf('Radiation pattern'), FontSize=15, ...
                    FontWeight='bold');
            end
            r = rad3D(k, d(n), phi, theta, dx(n,m));
            X = r.*cos(phi).*sin(theta);
            Y = r.*sin(phi).*sin(theta);
            Z = r.*cos(theta);
            surf(X, Y, Z);
            title(sprintf('d=%s, θ=%s', ddisp(n), angledisp(m)), ...
                FontSize=14);
            xlabel('x', FontSize=15, FontWeight='bold');
            ylabel('y', FontSize=15, FontWeight='bold');
            zlabel('z', FontSize=15, FontWeight='bold', Rotation=0);
            if axeq==2
                axis equal;
            end
        end
    end
end

% Συνάρτηση υπολογισμού του στερεού ακτινοβολίας σύμφωνα με την θεωρητική
% ανάλυση
function E = rad3D(k, d, phi, theta, dx)
    E = abs(cos(pi*cos(theta)/2)./sin(theta)).*abs(sin(8*k*d*cos(phi).*...
        sin(theta)+8*dx)./sin(k*d*cos(phi).*sin(theta)/2+dx/2)).* ...
            abs(sin(6*k*d*cos(theta))./sin(k*d*cos(theta)/2));
end

