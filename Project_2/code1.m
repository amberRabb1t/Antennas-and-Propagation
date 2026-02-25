% Κλείσιμο όλων των ανοιχτών figures και καθαρισμός του workspace, ώστε να
% σιγουρευτεί ότι καμία προηγούμενη εκτέλεση κώδικα δεν θα επηρεάσει την
% επόμενη. Επιπλέον, καθαρισμός του command window για ευκολότερη ανάγνωσή
% οποιουδήποτε μηνύματος τυπωθεί εκεί από το πρόγραμμα
close all;
clear;
clc;

% Ορισμός μήκους κύματος και κεντρικής συχνότητας 
% [Επώνυμο (Τ-Ω) --> λ = 2m]
c0 = 2.998*10^8;
wvl = 2;
f0 = (c0/wvl)/10^6;

% Ορισμός των μεταβλητών που δέχεται το πρόγραμμα: μήκος l διπόλων, ακτίνα
% r των wires, γωνία θ0, απόσταση d, αριθμός segments και αριθμός
% βημάτων για το frequency sweep
r = wvl/400;
l = wvl/2;
th0 = 30;
d = wvl/20;
segs = 10;
sweepstep = 255;

% Υπολογισμοί για τις συντεταγμένες των wires (μετατροπή από σφαιρικές σε
% καρτεσιανές για εισαγωγή στο αρχείο .nec)
phi = [0 60 120 180 240 300];
x = (l/2)*cosd(phi)*sind(th0);
y = (l/2)*sind(phi)*sind(th0);
z = d/2+(l/2)*cosd(th0);

% Δημιουργία του .nec αρχείου σύμφωνα με το συντακτικό του NEC
fileID = fopen('biconical.nec','w');

fprintf(fileID,'CM Biconical antenna\r\nCM Total cone opening angle %g deg.\r\nCE\r\n', 2*th0);

tag = 1;
for n=1:length(phi)
    fprintf(fileID,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\r\n', tag, segs, 0, 0, d/2, x(n), y(n), z, r);
    tag = tag + 1;
end

fprintf(fileID,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\r\n', tag, 1, 0, 0, -d/2, 0, 0, d/2, r);

for n=1:length(phi)
    tag = tag + 1;
    fprintf(fileID,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\r\n', tag, segs, 0, 0, -d/2, x(n), y(n), -z, r);
end

fprintf(fileID,'GE 0\r\n');
fprintf(fileID,'EX 0 %d 1 0 1 0\r\n',length(phi)+1);
fprintf(fileID,'FR 0 %d 0 0 %g %g\r\nEN\r\n',sweepstep,0.5*f0,(4*f0-0.5*f0)/sweepstep);

fclose(fileID);

