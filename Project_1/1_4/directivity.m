% Συνάρτηση υπολογισμού της κατευθυντικότητας μη-ομοιόμορφης
% στοιχειοκεραίας, σύμφωνα με την εκφώνηση και την θεωρητική ανάλυση
function D = directivity(p)
    p(6) = p(5);
    p(7) = p(4);
    p(8) = p(3);
    p(9) = p(2);
    p(10) = p(1);
    numerator = 0;
    for n=0:9
        numerator = numerator + p(n+1);
    end
    numerator = numerator^2;
    denominator = 0;
    for n=0:9
        for m=0:9
            denominator = denominator + p(n+1)*p(m+1)*sinc((n-m));
        end
    end
    D = -numerator/denominator;
end

