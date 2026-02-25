% Κλείσιμο όλων των ανοιχτών figures και καθαρισμός του workspace, ώστε να
% σιγουρευτεί ότι καμία προηγούμενη εκτέλεση κώδικα δεν θα επηρεάσει την
% επόμενη. Επιπλέον, καθαρισμός του command window για ευκολότερη ανάγνωσή
% οποιουδήποτε μηνύματος τυπωθεί εκεί από το πρόγραμμα
close all;
clear;
clc;

% Ορισμός ρυθμίσεων για τον γενετικό αλγόριθμο, σύμφωνα με την εκφώνηση
options = optimoptions('ga', PopulationSize=200, MaxGenerations=1000, ...
    MaxStallGenerations=500, MaxStallTime=200, PlotFcn='gaplotbestf');

[I,D]=ga(@directivity,5,[],[],[],[],[1 1 1 1 1],[4 4 4 4 4],[],[],options);

disp(I);
disp(-D);

