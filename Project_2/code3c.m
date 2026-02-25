% Κλείσιμο όλων των ανοιχτών figures και καθαρισμός του workspace, ώστε να
% σιγουρευτεί ότι καμία προηγούμενη εκτέλεση κώδικα δεν θα επηρεάσει την
% επόμενη. Επιπλέον, καθαρισμός του command window για ευκολότερη ανάγνωσή
% οποιουδήποτε μηνύματος τυπωθεί εκεί από το πρόγραμμα
close all;
clear;
clc;

% Ορισμός μήκους κύματος και κεντρικής συχνότητας 
% [Επώνυμο (Τ-Ω) --> f0 = 3GHz = 3*10^3MHz]
c0 = 2.998*10^8;
f0 = 3*10^3;
wvl = c0/(f0*10^6);
k = 2*pi/wvl;

% Ορισμός των μεταβλητών που δέχεται το πρόγραμμα: μήκος l διπόλων, ακτίνα
% r των wires, απόστασεις d και h μεταξύ τους και αριθμός segments
r = wvl/2000;
l = wvl/2;
d = 0.05*wvl:0.05*wvl:2*wvl;
h = 0.05*wvl:0.05*wvl:wvl;
segs = 21;

% Ορισμός πίνακα Zm για την αποθήκευση των τιμών αμοιβαίας σύνθετης 
% αντίστασης
Zm = zeros(length(h),length(d));

% Παράμετροι που χρειάζονται για την ανάκτηση των ζητούμενων πληροφοριών
% (ρεύματα στα δίπολα) από το εκάστοτε αρχείο .out (text file), που 
% δημιουργείται καλώντας τον nec solver για το εκάστοτε αρχείο .nec.
% Προκύπτουν από το formatting που καταλήγει να έχει το .out text file για
% την περίπτωση (όσον αφορά τον αριθμό segments) που έχουμε. 
% Ουσιαστικά, σε κάποιο section του αρχείου βρίσκονται οι πληροφορίες για 
% τα ρεύματα, formatted σαν πίνακας με 10 στήλες και γραμμές ίσες με το 
% διπλάσιο του αριθμού των segments (διότι έχουμε 2 δίπολα). 
% Άρα το formatSpec (%f) που θα χρησιμοποιηθεί στην έτοιμη συνάρτηση Matlab
% textscan() πρέπει να εφαρμοστεί 10*2*segs=20*segs φορές για να λάβουμε 
% όλον τον πίνακα.
% Ο "πίνακας" αυτός εμφανίζεται σε διαφορετικό line του text file ανάλογα 
% τα segments. Για 21 segments βρέθηκε ότι ξεκινά στο 154ο line
readreps = 20*segs;
lineskip = 153;

% Υπολογισμός της αυτοαντίστασης κεραίας λ/2 χρησιμοποιώντας τους
% θεωρητικούς τύπους
kl = k*l;
g = double(eulergamma);
a = wvl*10^(-5);
Rs = 60*(g+log(kl)-cosint(kl))+30*sin(kl)*(sinint(2*kl)-2*sinint(kl))+ ...
    30*cos(kl)*(g+log(kl/2)+cosint(2*kl)-2*cosint(kl));
Xs = 60*sinint(kl)+30*cos(kl)*(2*sinint(kl)-sinint(2*kl))-30*sin(kl)* ...
    (2*cosint(kl)-cosint(2*kl)-cosint(2*k*a^2/l));
Zs = Rs + 1i*Xs;

% Διπλή λούπα που, για όλες τις περιπτώσεις απόστασης h, διατρέχει τις
% περιπτώσεις αποστάσης d δημιουργώντας κάθε φορά το αντίστοιχο .nec file,
% για το οποίο καλεί τον nec solver. Στη συνέχεια διαβάζει τις πληροφορίες 
% από το αρχείο .out που δημιουργείται, με τρόπο που εξηγήθηκε παραπάνω. 
% Από τον "πίνακα" που αντλείται μας ενδιαφέρουν μόνο οι στήλες 7 και 8, 
% που έχουν τα ρεύματα. Οι τιμές αυτών αθροίζονται για να βρεθούν τα
% ρεύματα, από τα οποία υπολογίζεται και αποθηκεύεται η εκάστοτε Zm (οι 
% γραμμές είναι formatted έτσι ώστε οι πρώτες μισές να περιέχουν πληροφορία
% για το 1ο δίπολο και οι δεύτερες μισές για το 2ο -- έτσι βρίσκεται το I1 
% και το I2)
for v=1:length(h)
    for m=1:length(d)
        fileID = fopen('C:\4nec2\out\echelon.nec','w');

        fprintf(fileID,'CM Two in-echelon half-wavelength dipoles.\r\nCE\r\n');
        fprintf(fileID,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\r\n', 1, segs, d(m)/2, 0, 0, d(m)/2, 0, l, r);
        fprintf(fileID,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\r\n', 2, segs, -d(m)/2, 0, h(v), -d(m)/2, 0, h(v)+l, r);
        fprintf(fileID,'GE 0\r\n');
        fprintf(fileID,'EX 0 1 %d 0 1 0\r\n',int64(segs/2));
        fprintf(fileID,'FR 0 1 0 0 %g 0\r\n',f0);
        fprintf(fileID,'RP 0 37 73 1003 -180 0 5 5 0 0\r\nEN\r\n');

        fclose(fileID);

        system('(echo C:\4nec2\out\echelon.nec & echo C:\4nec2\out\echelon.out) | C:\4nec2\exe\nec2dxs1k5.exe');
        clc;

        fileID = fopen('C:\4nec2\out\echelon.out', 'r');
        currents = cell2mat(textscan(fileID,'%f', readreps, 'HeaderLines',lineskip));
        fclose(fileID);

        I1 = 0;
        I2 = 0;
        cnt = 0;
        for u=1:segs
            I1 = I1 + currents(7+cnt) + 1i*currents(8+cnt);
            cnt = cnt + 10;
        end
        for u=1:segs
            I2 = I2 + currents(7+cnt) + 1i*currents(8+cnt);
            cnt = cnt + 10;
        end

        Zm(v,m) = -Zs*(I2/I1);
    end
end

% Ορισμός των d και h ως συντεταγμένες ενός grid, ώστε να μπορέσουν να
% δημιουργηθούν τα 3D surface plots. Στη συνέχεια, δημιουργία των plots
[d,h] = meshgrid(d,h);
figure;
surf(d/wvl,h/wvl,real(Zm));
title('R_m vs d/λ & h/λ', FontSize=15);
zlabel('R_m (ohms)', FontSize=14, FontWeight='bold');
xlabel('Distance over wavelength d/λ', FontSize=14, FontWeight='bold');
ylabel('Distance over wavelength h/λ', FontSize=14, FontWeight='bold');
colorbar;
figure;
surf(d/wvl,h/wvl,imag(Zm));
title('X_m vs d/λ & h/λ', FontSize=15);
zlabel('X_m (ohms)', FontSize=14, FontWeight='bold');
xlabel('Distance over wavelength d/λ', FontSize=14, FontWeight='bold');
ylabel('Distance over wavelength h/λ', FontSize=14, FontWeight='bold');
colorbar;

