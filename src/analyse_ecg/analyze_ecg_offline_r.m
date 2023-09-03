%#eml
 function[rpeaks, art, correlation] = analyze_ecg_offline_r(signal, abtastrate)
 
 % notiz: IEM-abtastrate = 200
 
 if size(signal, 2) ~= 1
     error('Signal muss eindimensionaler Spaltenvektor sein');
 end

% Vektoren vorbereiten
differenz = zeros(1, length(signal));

origsignal = signal;
smoothsignal = smooth_emlc(signal, abtastrate/5);
signal = signal - smoothsignal(:,1); %!!!!!!!!!

% Fenstergröße zur Berechnung der Standardabweichung
standardabweichungsfenster = floor(0.2 * abtastrate);

% Lokale Standardabweichung zur Erkennung kurzfristiger Bewegungsartefakte
standardabweichung = movingstd(signal, standardabweichungsfenster)';

% Erkennung von Rauschen
% nulllinienfenster = floor(1.5 * abtastrate);
% for i = 1:length(signal)
%     a = max(1, i-nulllinienfenster);
%     o = min(length(signal), i+nulllinienfenster);
%     bereich = signal(a:o);
%     iqr_range(i) = iqr(bereich)/(max(bereich) - min(bereich));
%     nulldurchgang(i) = sum(abs(diff(sign(bereich - mean(bereich)))));
%     rauschen(i) = iqr_range(i) * nulldurchgang(i);
%     kurt(i) = kurtosis(bereich);
% end



% 2-Werte-Differenz (diskrete erste Ableitung) des Signals bilden
differenz(2:end-1) = signal(3:end) - signal(1:end-2);

% Fenstergröße für Amplitudenbestimmung festlegen
windowsize = floor(0.028 * abtastrate);
% windowsize = floor(0.05 * abtastrate);

% Amplituden von Signal und Differenzen berechnen
[runMin,runMax] = runningExtreme(signal, 2*windowsize+1, 'both');
signalampl = runMax - runMin;
[runMin,runMax] = runningExtreme(differenz, 2*windowsize+1, 'both');
diffampl = runMax - runMin;


% Featuresignal aus Amplituden von Signal und erster Ableitung
% zusammensetzen (Signal stärker gewichtet durch Quadrierung)
featuresignal = signalampl .* signalampl .* diffampl;


% Lokales Maximum des Featuresignals bilden (= eindimensionale Dilatation)
% -> Verbreiterung der Suchbereiche
maxwindow = floor(0.1 * abtastrate);
featuresignal = runningExtreme(featuresignal,2*maxwindow+1,'max');


% maxwindow = floor(0.15 * abtastrate);
% [runMin,runMax] = runningExtreme(signal .* signal .* differenz', 2*maxwindow+1, 'both');
% maxfeature = 4 * (runMax - runMin);
% % maxfeature = sqrt(maxfeature)*sqrt(max(maxfeature));

% maxfeature = 50*smooth_emlc(abs(signal .* signal .* differenz'), 2*windowsize+1);
% [runMin,runMax] = runningExtreme(maxfeature, 2*maxwindow+1, 'both');
% maxfeature = runMax - runMin;

% Variablen Grenzwert aus lokalem Mittelwert des Featuresignals und
% Standardabweichung des Signals berechnen
limit = smooth_emlc(featuresignal, 3*abtastrate);
limit = limit * 0.5 + 3*(standardabweichung').^3;

% Detektieren, wann das Featuresignal über dem Grenzwert liegt
feature = featuresignal > limit;


% Anfangs- und Endindizes der Suchbereiche bestimmen
feature = diff(feature);
start = find(feature == 1);
stop = find(feature == -1);


% Randbereiche überprüfen
if ~isempty(start) && ~isempty(stop)
    % Erster Startwert nach erstem Stopwert
    if start(1) > stop(1)
        start = [1; start];
    end
    
    % Letzer Stopwert vor letztem Startwert
    if stop(end) < start(end)
        stop = [stop; length(signal)];
    end
    
    % Führt sonst zu fehlerhaftem Matrixzugriff
    if start(1) <= 0
        start(1) = 1;
    end

    % Führt sonst zu fehlerhaftem Matrixzugriff
    if stop(end) > length(signal)
        stop(end) = length(signal);
    end
end


% Anzahl der Suchbereiche
number = min(length(start), length(stop));


% Init Ergebnisse nach dem ersten Durchlauf
firstpass = zeros(1, number);

% Für die Qualitätsbestimmung
testsignal_smooth = smooth_emlc(signal, (abtastrate/50)+1);

% R-Peak in jedem Suchbereich bestimmen
for j = 1:number
    % Werteverteilung im 3-Sekunden-Bereich um R-Zacke: Kurtosis
    a = max(1, start(j) - floor(2 * abtastrate));
    e = min(length(signal), stop(j) + floor(2 * abtastrate));
    
    kurt = kurtosis_emlc(signal(a:e), 0, 1);
    
    testsignal = testsignal_smooth(a:e);
    m = sum(testsignal)/length(testsignal);
    nulldurchgang = 2 * sum(abs(diff(testsignal > m))) / (e-a);
    
    s = sort(origsignal(a:e));
    resolution = sum(diff(s)~=0)+1;
    
    
%     o = origsignal(a:e);
%     d = abs(diff(o));
%     bla = range(o)/min(d(d>0));
    
    if (nulldurchgang < 0.3 || kurt > 4) && resolution > 30
        % Suchbereich vergrößern, falls zu klein
        if (stop(j) - start(j)) < 0.15 * abtastrate
            stop(j) = min(length(signal), floor(start(j) + 0.15 * abtastrate));
        else

        % Peak im Suchbereich finden        
        [ignore, ind] = findbigpeak(origsignal(start(j):stop(j)));
        firstpass(j) = ind + start(j) - 1;
        end
    end
end


% Bereiche ohne Peak entfernen
start = start(firstpass > 0);
stop = stop(firstpass > 0);
firstpass = firstpass(firstpass > 0);
number = length(firstpass);


% 80ms - Fenstergröße für den Suchbereich von QRS-Komplexen
qrslength = floor(0.080 * abtastrate);

% Array für QRS-Komplexe vorbereiten
parts = zeros(number, 2 * qrslength + 1);
extrasystolen = zeros(number, 2 * qrslength + 1);

% QRS-Komplexe separieren
for i = 1:number
    if firstpass(i)-qrslength > 0 && firstpass(i)+qrslength <= length(signal)
        parts(i,:) = signal(firstpass(i)-qrslength : firstpass(i)+qrslength);
        median_val = median(parts(i,:));
        median_val = median_val(1);
        parts(i,:) = parts(i,:) - median_val;
    end
end


% Prüfen, ob genug QRS-Komplexe gefunden wurden (wenn nicht, abbruch)
if size(parts, 1) < 2
    rpeaks = 0;
    featuresignal = 0;
    limit = 0;
    art = 0;
    correlation = 0;
    return
end

% Prüfen, ob die Auflösung des Signals groß genug ist
s = sort(origsignal);
resolution = sum(diff(s)~=0)+1;
if resolution < 70
    rpeaks = 0;
    featuresignal = 0;
    limit = 0;
    art = 0;
    correlation = 0;
    return
end

% Find 75% of the parts with higher positive peak
peakvalues = max(parts, [], 2);
rpeaks = zeros(2*number, 1);
correlation = nan(2*number, 1);
[v i] = sort(peakvalues);
half = ceil(length(i)/4);
selectedparts = parts(i(half:end),:);

% !! Template = Median der gefundenen QRS-Komplexe
% Find most representative beat as template
[v i] = max(sum(corrcoef(selectedparts')));
template = selectedparts(i,:); %median(parts(i(half:end),:), 1);
% template = median(parts, 1);

templateampl = template(ceil(length(template)/2)) - mean(template);

% In jedem Suchbereich Korrelation von Signal mit Template bilden:
% Ergebnispunkte verschieben, wenn in naher Umgebung eine bessere
% Korrelation auftritt bzw. entfernen, wenn die Korrelation generell zu
% klein ist.
lastrr = 0;
extrasystolenindex = zeros(number, 1);
rpeakcounter = 1;
extrasystolencounter = 1;

for i = 1:number
    % Korrelation im Suchbereich
    if (start(i)-qrslength > 0) && (stop(i)+qrslength < length(signal))
        
        korrelation = movingCorr(template, signal((start(i)-qrslength):(stop(i)+qrslength)));
        
% Version for coder:
%         for j = 1:(stop(i)-start(i)+1)
%             if j > maxsearchlength
%                 break
%             end
%             korrelation(j) = corr_emlc(template', signal_slide(:,j));
%         end
        
        % Peak in der Korrelation suchen
        [korrval, ind] = findsmallpeak(korrelation);
%         [korrval, ind] = findsmallpeak(combineddetection);
%         curamplitudefraction = amplitudefractions(ind);
%         curamplitudeerror = abs(1 - curamplitudefraction);
        ind = ind + start(i) - 1;

        % Wenn Korrelation groß genug: lokalen Peak im Signal als Ergebnis
        % speichern
        % Mindestens nötige Korrelation abhäng von der Veränderung des
        % RR-Intervalls
        if lastrr ~= 0 && (ind - rpeaks(rpeakcounter-1)) < 0.7*lastrr; %median(diff(rpeaks(rpeaks>0)))
            minkorr = 0.8;
        else
            minkorr = 0.7;
        end

        % Korrelation groß genug für normale R-Zacke und Herzfrequenz
        % kleiner als 250
        if korrval > minkorr && (lastrr == 0 || (ind - rpeaks(rpeakcounter-1)) > ((60/250) * abtastrate))% && curamplitudeerror < 0.7
            [ignore, fin] = findsmallpeak(origsignal(ind-8:ind+8));
            fin = fin + ind - 9;
            
            if fin > qrslength
                rpeaks(rpeakcounter) = fin;
                correlation(rpeakcounter) = korrval;
                if rpeakcounter > 1
                    lastrr = rpeaks(rpeakcounter) - rpeaks(rpeakcounter-1);
                end
                rpeakcounter = rpeakcounter + 1;
            end
            
        % Als mögliche Extrasystole vormerken
        else
            if extrasystolencounter <= number
                extrasystolen(extrasystolencounter, :) = parts(i,:);
                extrasystolenindex(extrasystolencounter) = firstpass(i);
                extrasystolencounter = extrasystolencounter + 1;
            end
        end
    end
end

extrasystemplate = zeros(1, 2 * qrslength + 1);

% Extrasystolen bestimmen
if extrasystolencounter > 2
    extrasystolen = extrasystolen(1:extrasystolencounter-1, :);
    % Find most representative beat as template
    [v i] = max(sum(corrcoef(extrasystolen')));
    extrasystemplate = extrasystolen(i,:); %median(parts(i(half:end),:), 1);
%     extrasystemplate = median(extrasystolen);
    
    for i = 1 : extrasystolencounter-1
        korrval = corr_emlc(extrasystolen(i,:)', extrasystemplate');
        if korrval < 0.6
            extrasystolenindex(i) = 0;
        end
    end
else
    extrasystolenindex(1) = 0;
end


% Aufhören, wenn keine neuen R-Zacken mehr gefunden werden
maxlen = length(rpeaks);
repeat = 1;
notnullpeaks = zeros(1, 2*number);
while repeat == 1
    repeat = 0;
    
    % RR-Intervalle berechnen
    if (sum(rpeaks > 0) + sum(extrasystolenindex > 0)) <= 2*number
        notnullpeaks = [rpeaks(rpeaks > 0)', extrasystolenindex(extrasystolenindex > 0)'];
    end
    
    notnullpeaks = unique(notnullpeaks);
    
    if length(notnullpeaks) < 2
        break
    end
    
	rr_int = notnullpeaks(2:end) - notnullpeaks(1:end-1);
    
    if isempty(rr_int) %damit es nicht zu einem fehler kommt, wenn nur ein element in notnullpeaks ist.
        break
    end
    
    % Wenn Abstand zwischen zwei RR-Intervallen zu groß: Davon ausgehen, dass
    % eine R-Zacke übersehen wurde
    missing = rr_int(2:end) > 1.8 * rr_int(1:end-1) & rr_int(2:end) > ((60/125) * abtastrate);

    % Bereiche mit fehlender R-Zacke mit geringerer Mindestkorrelation absuchen
    missing_ind = find(missing == 1);
    for z = 1:length(missing_ind)
        i = missing_ind(z);
        
        index = ceil(2 * notnullpeaks(i+1) - notnullpeaks(i));

        if index-3*qrslength > 0 && index+3*qrslength <= length(signal)
            
            % Korrelation im Suchbereich
            korrelation = movingCorr(template, signal(index-3*qrslength:index+3*qrslength));

            % realive Amplitude im Vergleich zu Template
            windowmean = smooth_emlc(signal(index-3*qrslength : index+3*qrslength), 2*qrslength+1);
            windowmiddle = signal(index-2*qrslength : index+2*qrslength);
            amplitudefractions = (windowmiddle - windowmean(qrslength+1:end-qrslength))/templateampl;
        

            % Peak in der Korrelation suchen
            [korrval, ind] = findsmallpeak(korrelation);
            curamplitudefraction = amplitudefractions(ind);
            curamplitudeerror = abs(1 - curamplitudefraction);
            ind = ind + index - 2*qrslength - 1;
            
            if korrval > 0.6 && curamplitudeerror < 0.80  %&& featuresignal(ind) > limit(ind)/2
                [ignore, fin] = findsmallpeak(origsignal(ind-10:ind+10));
                
                if fin > 0 && rpeakcounter < maxlen
                    fin = fin + ind - 11;
                    
                    nearestrr = min(abs(rpeaks - fin));
                    
                    if (rpeakcounter + 1) < maxlen
                        if nearestrr > ((60/250) * abtastrate)
                            rpeaks(rpeakcounter) = fin;
                            correlation(rpeakcounter) = korrval;
                            rpeakcounter = rpeakcounter + 1;
                            repeat = 1;
                        end
                    else
                        break;
                    end
                end
            end
        end
    end
end


% R-Zacken neu zuordnen (eventuell sind ja Extrasystolen dabei)
if sum(extrasystolenindex > 0) > 1
    for i = 1:(rpeakcounter-1)
        p = signal(rpeaks(i)-qrslength:rpeaks(i)+qrslength);   
        ct = corr_emlc(p, template');
        ce = corr_emlc(p, extrasystemplate');
        if (ct < ce) && (ce > 0.7)% && (ct < 0.7)
            extrasystolenindex(extrasystolencounter) = rpeaks(i);
            extrasystolencounter = extrasystolencounter + 1;
            rpeaks(i) = 0;
        end
    end
end


correlation = correlation(rpeaks > 0);

rpeaks = rpeaks(rpeaks > 0);
[rpeaks, i] = unique(rpeaks');

correlation = correlation(i)';

extrasystolenindex = extrasystolenindex(extrasystolenindex > 0);
extrasystolenindex = unique(extrasystolenindex');

art = [zeros(1, length(rpeaks)), ones(1, length(extrasystolenindex))] + 1;

correlation = [correlation, nan(1, length(extrasystolenindex))];

rpeaks = [rpeaks, extrasystolenindex];
[rpeaks, i] = unique(rpeaks, 'first');
art = art(i);

correlation = correlation(i);

if sum(art) > 1.5 * length(art)
    art = mod(art, 2) + 1;
end

return


% Findet den Peak mit maximaler Amplitude
function[value, index] = findsmallpeak(data)
d = diff(data);
i = find(d(1:end-1) > 0 & d(2:end) <= 0) + 1;

if isempty(i)
%     value = -Inf; %data(1);
%     index = -Inf; %1;
    [value, index] = max(data);
else
    [value, maxind] = max(data(i));
    index = i(maxind);
end

% alter code, gleiche ergebnisse
% value = data(1);
% index = 1;
% for i = 2:length(data)-1
%     if data(i-1) < data(i) && data(i+1) <= data(i) && value < data(i)
%         index = i;
%         value = data(i);
%     end
% end

return


% Findet den Peak mit maximalem Berg-Tal-Verhältnis
function[value, index] = findbigpeak(data)

% Finde echte peaks -> Nulldurchgang in 1. Ableitung
t = diff(data);
s = sign(t);
s(s<0) = 0;
x = diff(s);
rising = find(x>0)+1;
falling = find(x<0)+1;

% Nichts gefunden
if isempty(falling)
    value = -Inf;
    index = -Inf;
    return
end

% Jeder positive Peak sollte 2 negative Peaks als Nachbarn haben
if isempty(rising) || rising(1) > falling(1)
    rising = [1; rising];
end

if length(rising) == length(falling)
    rising = [rising; length(data)];
end

% Positiven Peak mit größtem Berg-Tal-Verhältnis finden (fallende Kante
% etwas stärker gewichtet als steigende Kante -> bessere Erkennung von
% rS - Komplexen)
risingElevation = data(falling) - data(rising(1:end-1));
fallingElevation = data(falling) - data(rising(2:end));
elevation = risingElevation + 1.9 * fallingElevation; %1.5
[value, i] = max(elevation);
index = falling(i);

% Passenden Datenwert zum gedundenen Index holen
if index > 0
    value = data(index);
end

return