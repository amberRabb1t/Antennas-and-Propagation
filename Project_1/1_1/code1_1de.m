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

% Ορισμός των τιμών που εξετάζουμε, αριθμητικά για τις πράξεις και σε
% string για την σωστή εμφάνισή τους με σύμβολα (και όχι νούμερα) κατά την
% τύπωση στο command window
deltax = [0 -pi/2 -sqrt(3)*pi/2 -pi];
ndx = length(deltax);
angledisp = ["0°" "30°" "60°" "90°"];
dx = wvl / 2;

% δ) Λούπες που υλοποιούν με διπλό άθροισμα Riemann το διπλό ολοκλήρωμα της
% θεωρητικής ανάλυσης και τυπώνουν τα αποτελέσματα. Το βήμα είναι μισή
% μοίρα και τα φ,θ αναγκαστικά δεν ξεκινούν από το 0, ώστε να αποφευχθούν 
% απροσδιοριστίες
Rsum = zeros(1,ndx);
D = zeros(1,ndx);
fprintf('δ) Riemann sum method:\n\t');
for n=1:ndx
    for phi = 0.0000001:pi/360:2*pi
        for theta = 0.0000001:pi/360:pi
            Rsum(n) = Rsum(n) + calcE(k, dx, phi, theta, deltax(n))^2* ...
                sin(theta)*(pi/360)^2;
        end
    end
    D(n) = 36864*4*pi/Rsum(n);
    fprintf('θ=%s\t\t', angledisp(n));
end
fprintf('\n\t');
for n=1:ndx
    fprintf('D=%g\t', D(n));
end

% ε) Σχεδιασμός του στερεού ακτινοβολίας και υπολογισμός κατευθυντικότητας,
% έχοντας σχεδιάσει την κεραία να λειτουργεί ως ακροπυροδοτική
% Hansen-Woodyard. Ορίζεται το εύρος γωνιών σε μορφή συντεταγμένων, καθώς
% και οι παράμετροι dx=λ/4 και δx=-π/2-0.1825 όπως βρέθηκε στην θεωρητική
% ανάλυση. Η κατευθυντικότητα υπολογίζεται (υπολογιστικά) όμοια με πριν
for n=1:2
    [phi, theta] = meshgrid(0:0.009:2*pi, 0:0.009:pi);
    dx = wvl / 4;
    dz = wvl / 2;
    deltax = -pi/2 -0.1825;
    r = rad3DHW(k, dx, dz, phi, theta, deltax);
    X = r.*cos(phi).*sin(theta);
    Y = r.*sin(phi).*sin(theta);
    Z = r.*cos(theta);
    figure;
    surf(X, Y, Z);
    if n==2
         title(sprintf(['Radiation pattern — axis equal\nOperating as '...
             'Hansen-Woodyard end-fire']), FontSize=14);
    else
         title(sprintf(['Radiation pattern\nOperating as ' ...
             'Hansen-Woodyard end-fire']), FontSize=14);
    end
    xlabel('x', FontSize=15, FontWeight='bold');
    ylabel('y', FontSize=15, FontWeight='bold');
    zlabel('z', FontSize=15, FontWeight='bold', Rotation=0);
    if n==2
        axis equal;
    end
end

RsumH = 0;
fprintf(['\n\nε) Operating as Hansen-Woodyard end-fire\n\tRiemann sum' ...
    ' method:\n']);
for phi = 0.0000001:pi/360:2*pi
    for theta = 0.0000001:pi/360:pi
        RsumH = RsumH +calcEHW(k, dx, dz, phi, theta, deltax)^2* ...
            sin(theta)*(pi/360)^2;
    end
end
DH = 17129.993*4*pi/RsumH;
fprintf('\tD=%g\n', DH);


% Συνάρτηση υπολογισμού του μέτρου του Ε (χωρίς Ε0) σύμφωνα με την
% θεωρητική ανάλυση
function E = calcE(k, d, phi, theta, dx)
    E = abs(cos(pi*cos(theta)/2)/sin(theta))*abs(sin(8*k*d*cos(phi)* ...
            sin(theta)+8*dx)/sin(k*d*cos(phi)*sin(theta)/2+dx/2))* ...
                abs(sin(6*k*d*cos(theta))/sin(k*d*cos(theta)/2));
end

% Συνάρτηση υπολογισμού του μέτρου του Ε (χωρίς Ε0) σύμφωνα με την
% θεωρητική ανάλυση (Hansen-Woodyard)
function E = calcEHW(k, dx, dz, phi, theta, deltax)
    E = abs(cos(pi*cos(theta)/2)/sin(theta))*abs(sin(8*k*dx*cos(phi)* ...
            sin(theta)+8*deltax)/sin(k*dx*cos(phi)*sin(theta)/2+deltax ...
            /2))*abs(sin(6*k*dz*cos(theta))/sin(k*dz*cos(theta)/2));
end

% Συνάρτηση υπολογισμού του στερεού ακτινοβολίας σύμφωνα με την θεωρητική
% ανάλυση (Hansen-Woodyard)
function E = rad3DHW(k, dx, dz, phi, theta, deltax)
    E = abs(cos(pi*cos(theta)/2)./sin(theta)).*abs(sin(8*k*dx*cos(phi).*...
        sin(theta)+8*deltax)./sin(k*dx*cos(phi).*sin(theta)/2+deltax/2 ...
        )).*abs(sin(6*k*dz*cos(theta))./sin(k*dz*cos(theta)/2));
end

