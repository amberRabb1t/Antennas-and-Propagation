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
% r των wires, απόσταση d μεταξύ τους και αριθμός segments
r = wvl/2000;
l = wvl/2;
d = 0.05*wvl:0.05*wvl:3*wvl;
segs = [5 21];

% Ορισμός πίνακα Zm για την αποθήκευση των τιμών αμοιβαίας σύνθετης 
% αντίστασης
Zm = zeros(length(d),length(segs));

% Παράμετροι που χρειάζονται για την ανάκτηση των ζητούμενων πληροφοριών
% (ρεύματα στα δίπολα) από το εκάστοτε αρχείο .out (text file), που 
% δημιουργείται καλώντας τον nec solver για το εκάστοτε αρχείο .nec.
% Προκύπτουν από το formatting που καταλήγει να έχει το .out text file για
% τις δύο περιπτώσεις που έχουμε. 
% Ουσιαστικά, σε κάποιο section του αρχείου βρίσκονται οι πληροφορίες για 
% τα ρεύματα, formatted σαν πίνακας με 10 στήλες και γραμμές ίσες με το 
% διπλάσιο του αριθμού των segments (διότι έχουμε 2 δίπολα). 
% Άρα το formatSpec (%f) που θα χρησιμοποιηθεί στην έτοιμη συνάρτηση Matlab
% textscan() πρέπει να εφαρμοστεί 10*2*segs=20*segs φορές για να λάβουμε 
% όλον τον πίνακα.
% Ο "πίνακας" αυτός εμφανίζεται σε διαφορετικό line του text file ανάλογα 
% τα segments. Για 5 segments βρέθηκε ότι ξεκινά στο 122ο line ενώ για 21
% segments ξεκινά στο 154ο line
readreps = 20*segs;
lineskip = [121 153];

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

% Διπλή λούπα που, για τις 2 περιπτώσεις αριθμού segments, διατρέχει τις
% περιπτώσεις αποστάσης d δημιουργώντας κάθε φορά το αντίστοιχο .nec file,
% για το οποίο καλεί τον nec solver. Στη συνέχεια διαβάζει τις πληροφορίες 
% από το αρχείο .out που δημιουργείται, με τρόπο που εξηγήθηκε παραπάνω. 
% Από τον "πίνακα" που αντλείται μας ενδιαφέρουν μόνο οι στήλες 7 και 8, 
% που έχουν τα ρεύματα. Οι τιμές αυτών αθροίζονται για να βρεθούν τα
% ρεύματα, από τα οποία υπολογίζεται και αποθηκεύεται η εκάστοτε Zm (οι 
% γραμμές είναι formatted έτσι ώστε οι πρώτες μισές να περιέχουν πληροφορία
% για το 1ο δίπολο και οι δεύτερες μισές για το 2ο -- έτσι βρίσκεται το I1 
% και το I2)
for n=1:length(segs)
    for m=1:length(d)
        fileID = fopen('C:\4nec2\out\parallel.nec','w');

        fprintf(fileID,'CM Two parallel half-wavelength dipoles.\r\nCE\r\n');
        fprintf(fileID,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\r\n', 1, segs(n), d(m)/2, 0, l/2, d(m)/2, 0, -l/2, r);
        fprintf(fileID,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\r\n', 2, segs(n), -d(m)/2, 0, l/2, -d(m)/2, 0, -l/2, r);
        fprintf(fileID,'GE 0\r\n');
        fprintf(fileID,'EX 0 1 %d 0 1 0\r\n',int64(segs(n)/2));
        fprintf(fileID,'FR 0 1 0 0 %g 0\r\n',f0);
        fprintf(fileID,'RP 0 37 73 1003 -180 0 5 5 0 0\r\nEN\r\n');

        fclose(fileID);

        system('(echo C:\4nec2\out\parallel.nec & echo C:\4nec2\out\parallel.out) | C:\4nec2\exe\nec2dxs1k5.exe');
        clc;

        fileID = fopen('C:\4nec2\out\parallel.out', 'r');
        currents = cell2mat(textscan(fileID,'%f', readreps(n), 'HeaderLines',lineskip(n)));
        fclose(fileID);

        I1 = 0;
        I2 = 0;
        cnt = 0;
        for u=1:segs(n)
            I1 = I1 + currents(7+cnt) + 1i*currents(8+cnt);
            cnt = cnt + 10;
        end
        for u=1:segs(n)
            I2 = I2 + currents(7+cnt) + 1i*currents(8+cnt);
            cnt = cnt + 10;
        end

        Zm(m,n) = -Zs*(I2/I1);
    end
end

% Δημιουργία των ζητούμενων διαγραμμάτων αμοιβαίας σύνθετης αντίστασης, για
% τις 2 περιπτώσεις αριθμού segments
for u=1:length(segs)
    figure(u);
    plot(d/wvl,real(Zm(:,u)));
    hold on;
    plot(d/wvl,imag(Zm(:,u)));
    yline(0,'--');
    title('Z_m vs d/λ', FontSize=15);
    ylabel('Mutual impedance Z_m (ohms)', FontSize=14, FontWeight='bold');
    xlabel('Distance over wavelength d/λ', FontSize=14, FontWeight='bold');
    legend('R_m','X_m');
end

