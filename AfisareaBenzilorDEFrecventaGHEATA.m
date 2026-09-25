%% Afisarea benzilor de frecventa
clear; clc; close all;
[fisier, cale] = uigetfile('*.set', 'Selecteaza fisierul EEG');
if isequal(fisier, 0)
    disp('Operatiune anulata de utilizator.');
    return;
end
EEG = pop_loadset('filename', fisier, 'filepath', cale);
data = EEG.data(1, :); % Procesăm primul canal 

% 2.PARAMETRI DE BAZĂ
Fs = EEG.srate;            % Frecvența de eșantionare
L = length(data);          % Lungimea semnalului
t = (0:L-1)*(1/Fs);        % Vectorul de timp
timp_afisare = 5;          % Secunde pentru vizualizarea formei de undă
NFFT = 2^nextpow2(L);      % Lungime pentru FFT (putere a lui 2)
f = Fs/2 * linspace(0, 1, NFFT/2); % Vectorul de frecvență

% Setări filtre 
Dpass = 0.0575; 
Dstop = 0.0001;
dens  = 20;
clear; clc; close all;

[fisier, cale] = uigetfile('*.set', 'Selecteaza fisierul EEG');
if isequal(fisier, 0)
    disp('eroare');
    return;
end

EEG = pop_loadset('filename', fisier, 'filepath', cale);
data = EEG.data(1, :); 

Fs = EEG.srate;        % Frecvența de eșantionare preluată din aparat
L = length(data);      % Lungimea semnalului
t = (0:L-1)*(1/Fs);    % Vectorul de timp
timp_afisare = 5;      % Fereastra de vizualizare (secunde)

NFFT = 2^nextpow2(L);  
f = Fs/2*linspace(0, 1, NFFT/2); 

% BANDA DELTA (0.5 - 4 Hz)
Fpass = 4;              
Fstop = 4.5;              
Dpass = 0.057501127785; 
Dstop = 0.0001;         
dens  = 20;             

[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
b1 = firpm(N, Fo, Ao, W, {dens});
Hd1 = dfilt.dffir(b1);
x1 = filter(Hd1, data);

h4 = figure; plot(t, x1, 'r'); 
title('Forma de undă pentru banda DELTA');
xlabel('Timp (secunde)'); ylabel('Amplitudine (\muV)');
xlim([0 timp_afisare]); 

Y1 = fft(x1, NFFT)/L;

h8 = figure; 
plot(f, 2*abs(Y1(1:NFFT/2)), 'r'); 
title('Spectrul de amplitudine unilateral pentru banda DELTA');
xlabel('Frecvență (Hz)'); ylabel('|Amplitudine|');
xlim([0 40]); 

% BANDA THETA (4-7 Hz)
Fstop1 = 3.5;          
Fpass1 = 4;             
Fpass2 = 7;             
Fstop2 = 7.5;           
Dstop1 = 0.0001;        
Dpass  = 0.057501127785;
Dstop2 = 0.0001;        

[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
b2 = firpm(N, Fo, Ao, W, {dens});
Hd2 = dfilt.dffir(b2);
x2 = filter(Hd2, data);

h5 = figure; plot(t, x2, 'm'); 
title('Forma de undă pentru banda THETA');
xlabel('Timp (secunde)'); ylabel('Amplitudine (\muV)');
xlim([0 timp_afisare]); 

Y2 = fft(x2, NFFT)/L;

h9 = figure; 
plot(f, 2*abs(Y2(1:NFFT/2)), 'm'); 
title('Spectrul de amplitudine unilateral pentru banda THETA');
xlabel('Frecvență (Hz)'); ylabel('|Amplitudine|');
xlim([0 40]); 

% BANDA ALPHA (8-12 Hz)
Fstop1 = 7.5;           
Fpass1 = 8;             
Fpass2 = 12;            
Fstop2 = 12.5;          

[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
b3 = firpm(N, Fo, Ao, W, {dens});
Hd3 = dfilt.dffir(b3);
x3 = filter(Hd3, data);

h6 = figure; plot(t, x3, 'g'); 
title('Forma de undă pentru banda ALPHA');
xlabel('Timp (secunde)'); ylabel('Amplitudine (\muV)');
xlim([0 timp_afisare]); 

Y3 = fft(x3, NFFT)/L;

h10 = figure; 
plot(f, 2*abs(Y3(1:NFFT/2)), 'g'); 
title('Spectrul de amplitudine unilateral pentru banda ALPHA');
xlabel('Frecvență (Hz)'); ylabel('|Amplitudine|');
xlim([0 40]); 

% BANDA BETA (12-30 Hz)
Fstop1 = 11.5;          
Fpass1 = 12;            
Fpass2 = 30;            
Fstop2 = 30.5;          

[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
b4 = firpm(N, Fo, Ao, W, {dens});
Hd4 = dfilt.dffir(b4);
x4 = filter(Hd4, data);

h7 = figure; plot(t, x4, 'b'); 
title('Forma de undă pentru banda BETA');
xlabel('Timp (secunde)'); ylabel('Amplitudine (\muV)');
xlim([0 timp_afisare]); 

Y4 = fft(x4, NFFT)/L;

h11 = figure; 
plot(f, 2*abs(Y4(1:NFFT/2)), 'b'); 
title('Spectrul de amplitudine unilateral pentru banda BETA');
xlabel('Frecvență (Hz)'); ylabel('|Amplitudine|');
xlim([0 40]);