    %procesarea unui singur set de date 
    data = readtable('C:\Users\Bianca\Documents\CercetareIRA\GHEATA\BBT-E12-AAB038-2025-03-10_10-51-56\EEG.csv');
    eeg_data = table2array(data(:, 5:16))';
    disp(size(eeg_data)) 
    
    %%  Procesare persoane
cd = 'C:\Users\Bianca\Documents\CercetareIRA\GHEATA';
datasetList = dir(fullfile(cd, 'BBT*')); 
datasetList = datasetList([datasetList.isdir]); 
    
    for i = 1:length(datasetList)
        dataName = datasetList(i).name;
        fprintf('--- Procesăm subiectul %d: %s ---\n', i, dataName);
        
        % 2. Importul datelor (.CSV)
        cale_csv = fullfile(datasetList(i).folder, dataName, 'EEG.csv');
        temp_table = readtable(cale_csv);
        eeg_data = temp_table{:, 5:16}'; 
        
        EEG = pop_importdata('dataformat','array','nbchan',12,'data','eeg_data','srate',256);
        EEG.setname = dataName;
    
        % 3. FILTRU TRECE-SUS (High-Pass) la 1 Hz
        EEG = pop_eegfiltnew(EEG, 'locutoff', 1);
    
        % 4. FILTRU TRECE-JOS (Low-Pass) la 40 Hz  
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 40);
    
        % 5. FILTRU NOTCH la 50 Hz 
        EEG = pop_eegfiltnew(EEG, 'locutoff', 48, 'hicutoff', 52, 'revfilt', 1);
    
        % 6. CURĂȚARE SPIKE 
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 5, 'ChannelCriterion', 0.8, 'BurstCriterion', 20, 'WindowCriterion', 'off');
    
        % 7. RE-REFERENȚIERE (Average Reference)
        EEG = pop_reref(EEG, []);
    
        % 8. SALVARE AUTOMATĂ
        nume_salvare = [dataName 'Procesat.set'];
        EEG = pop_saveset(EEG, 'filename', nume_salvare, 'filepath', cd);
        ALLEEG = []; EEG = [];
    end
    
%% Analiza Benzilor de Frecvență -SET GHEATA

cale_folder = 'C:\Users\Bianca\Documents\CercetareIRA\GHEATA\'; 
lista_fisiere = dir(fullfile(cale_folder, '*Procesat.set'));
nume_electrozi = {'Fp1', 'Fp2', 'F7', 'F3', 'F4', 'F8', 'P7', 'P3', 'P4', 'P8', 'O1', 'O2'};
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

rezultate_matrice = zeros(length(lista_fisiere), 5);
nume_subiecti = cell(length(lista_fisiere), 1);

fprintf('Procesare GHEATA pentru %d subiecți...\n', length(lista_fisiere));

for i = 1:length(lista_fisiere)
    EEG = pop_loadset('filename', lista_fisiere(i).name, 'filepath', cale_folder);
  
    for ch = 1:EEG.nbchan
        EEG.chanlocs(ch).labels = nume_electrozi{ch};
    end
    
    % 3. Căutare coordonate și SALVARE pe disc 
    EEG = pop_chanedit(EEG, 'lookup', which('standard-10-5-cap385.elp'));
    EEG = pop_saveset(EEG, 'savemode', 'resave');
   
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, i);
  
    % 5. Analiza Spectrală
    [spectra, freqs] = spectopo(EEG.data, EEG.pnts, EEG.srate, 'plot', 'off');
    
    deltaIdx = find(freqs >= 1 & freqs <= 4);
    thetaIdx = find(freqs > 4 & freqs <= 8);
    alphaIdx = find(freqs > 8 & freqs <= 13);
    betaIdx  = find(freqs > 13 & freqs <= 30);
    gammaIdx = find(freqs > 30 & freqs <= 40);
    
    rezultate_matrice(i,1) = mean(10.^(mean(spectra(:, deltaIdx), 2)/10));
    rezultate_matrice(i,2) = mean(10.^(mean(spectra(:, thetaIdx), 2)/10));
    rezultate_matrice(i,3) = mean(10.^(mean(spectra(:, alphaIdx), 2)/10));
    rezultate_matrice(i,4) = mean(10.^(mean(spectra(:, betaIdx),  2)/10));
    rezultate_matrice(i,5) = mean(10.^(mean(spectra(:, gammaIdx), 2)/10));
    
    nume_subiecti{i} = lista_fisiere(i).name;
    fprintf('Subiect %d: %s procesat și stocat în ALLEEG.\n', i, lista_fisiere(i).name);
end
TabelGheata = array2table(rezultate_matrice, 'VariableNames', {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'});
TabelGheata.Fisier = nume_subiecti;
TabelGheata = TabelGheata(:, [6 1 2 3 4 5]); 
writetable(TabelGheata, fullfile(cale_folder, 'REZULTATE_FINALE_GHEATA.csv'));
eeglab redraw; 
%% Afisarea benzilor de frecventa
clear; clc; close all;
[fisier, cale] = uigetfile('*.set', 'Selecteaza fisierul EEG');
if isequal(fisier, 0)
    disp('Operatiune anulata de utilizator.');
    return;
end
EEG = pop_loadset('filename', fisier, 'filepath', cale);
data = EEG.data(1, :); % Procesăm primul canal 

% 2. PARAMETRI DE BAZĂ
Fs = EEG.srate;            % Frecvența de eșantionare
L = length(data);          % Lungimea semnalului
t = (0:L-1)*(1/Fs);        % Vectorul de timp
timp_afisare = 5;          % Secunde pentru vizualizarea formei de undă
NFFT = 2^nextpow2(L);      % Lungime pentru FFT (putere a lui 2)
f = Fs/2 * linspace(0, 1, NFFT/2); % Vectorul de frecvență

% Setări generale filtre 
Dpass = 0.0575; 
Dstop = 0.0001;
dens  = 20;

% DELTA BAND (0.5 - 4 Hz)
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

NFFT = 2^nextpow2(L);  % Optimizare pentru FFT
f = Fs/2*linspace(0, 1, NFFT/2); % Vectorul de frecvențe

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
%% cod entropie 
%  PhysioNet: https://physionet.org/content/sampen/1.0.0/
cale_folder = 'C:\Users\Bianca\Documents\CercetareIRA\GHEATA\';
lista_fisiere = dir(fullfile(cale_folder, '*.set'));

% Parametrii SampEn - valori standard din literatura
m = 2;   
r = 0.2;  
nume_electrozi = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','Fz','Pz'};

% Matricea rezultate: un SampEn per canal per subiect
rezultate_sampen = zeros(length(lista_fisiere), 12);
nume_subiecti = cell(length(lista_fisiere), 1);
for i = 1:length(lista_fisiere)
    EEG = pop_loadset('filename', lista_fisiere(i).name, 'filepath', lista_fisiere(i).folder);
    nume_subiecti{i} = lista_fisiere(i).name;

    for ch = 1:EEG.nbchan
        EEG.chanlocs(ch).labels = nume_electrozi{ch};
    end

    fprintf('\nCalculez SampEn pentru: %s\n', lista_fisiere(i).name);
    % Calcul SampEn pentru fiecare canal
    for ch = 1:EEG.nbchan
        semnal = double(EEG.data(ch, 1:min(3000, EEG.pnts)));
        e = sampen(semnal, m, r, 'chebychev');
        rezultate_sampen(i, ch) = e;
        
    end

    disp(['Subiect ', num2str(i), ' din ', num2str(length(lista_fisiere)), ' finalizat.']);
end
% Salvare CSV 
TabelSampEn = array2table(rezultate_sampen, 'VariableNames', nume_electrozi);
TabelSampEn.Fisier = nume_subiecti;
TabelSampEn = TabelSampEn(:, [13 1:12]);
writetable(TabelSampEn, fullfile(cale_folder, 'REZULTATE_SAMPEN_GHEATA.csv'));
fprintf('\n=== GATA! SampEn salvat in REZULTATE_SAMPEN_GHEATA.csv ===\n');
disp(TabelSampEn);
%% CLUSTERING NESUPERVIZAT - K-MEANS EEG GHEATA
% Sursa 1: MathWorks. Cluster Analysis Example
%          https://www.mathworks.com/help/stats/cluster-analysis-example.html
%Sursa 2:AI

% 1. INCARCARE DATE
tabel_date = readtable('C:\Users\Bianca\Documents\CercetareIRA\GHEATA\REZULTATE_SAMPEN_GHEATA.csv');
canale_eeg = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','Fz','Pz'};
X = tabel_date{:, canale_eeg};
fprintf('Date incarcate: %d subiecti x %d canale\n', size(X,1), size(X,2));

% 2. NORMALIZARE
X_norm = normalize(X);

% 3. PCA - reducere dimensionalitate pentru vizualizare 2D
[coeff, score, ~, ~, explained] = pca(X_norm);
fprintf('PC1: %.1f%% | PC2: %.1f%% | Total: %.1f%%\n', ...
        explained(1), explained(2), explained(1)+explained(2));

% 4. ELBOW METHOD - numar optim de clustere
figure('Name','Elbow Method','Position',[100 100 500 400]);
sumd_vec = zeros(1,6);
for k = 1:6
    [~,~,sumd] = kmeans(X_norm, k, 'Replicates', 10, 'Display', 'off');
    sumd_vec(k) = sum(sumd);
end
plot(1:6, sumd_vec, '-o', 'LineWidth', 2, 'MarkerFaceColor','b', 'MarkerSize', 8);
xlabel('Numar clustere (k)');
ylabel('Suma distante intra-cluster');
title('Elbow Method - Numar optim clustere');
xline(2,'--r','k=2','LineWidth',1.5);
grid on;

% 5. K-MEANS k=2 - Sanatos vs Bolnav
rng(42);
[idx_clustere, centroizi, sumdist] = kmeans(X_norm, 2, ...
    'Distance',   'sqeuclidean', ...
    'Replicates', 10, ...
    'Display',    'off');
fprintf('\n=== REZULTATE CLUSTERING ===\n');
fprintf('Cluster 1: %d subiecti\n', sum(idx_clustere==1));
fprintf('Cluster 2: %d subiecti\n', sum(idx_clustere==2));

% 6. SILHOUETTE SCORE
s = silhouette(X_norm, idx_clustere);
fprintf('Silhouette Score mediu: %.4f\n', mean(s));

% 7. PURITATE CLUSTERE - cat de bine separa sanatosi vs bolnavi
media_c1 = mean(X(idx_clustere==1, :));
media_c2 = mean(X(idx_clustere==2, :));
fprintf('\nMedia SampEn per cluster:\n');
fprintf('%-6s | %-12s | %-12s | %s\n', 'Canal','Cluster 1','Cluster 2','Diferenta');
fprintf('%s\n', repmat('-',1,45));
for j = 1:length(canale_eeg)
    fprintf('%-6s | %-12.4f | %-12.4f | %.4f\n', ...
            canale_eeg{j}, media_c1(j), media_c2(j), media_c1(j)-media_c2(j));
end

% 8. GRAFICE
figure('Name','Clustering Nesupervizat - GHEATA', ...
       'Position',[50 50 1400 500]);

% Subplot 1 - Clustere pe PCA
subplot(1,3,1);
gscatter(score(:,1), score(:,2), idx_clustere, ...
         [0.2 0.5 0.8; 0.8 0.3 0.2], 'os', 8);
hold on;
centroizi_pca = centroizi * coeff(:,1:2);
plot(centroizi_pca(:,1), centroizi_pca(:,2), 'kx', ...
     'MarkerSize',15, 'LineWidth',3);
hold off;
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('Clustere K-Means pe PCA');
legend({'Cluster 1','Cluster 2','Centroizi'}, 'Location','best','FontSize',8);
grid on;

% Subplot 2 - Media SampEn per canal per cluster
subplot(1,3,2);
bar([media_c1; media_c2]', 'grouped');
set(gca,'XTickLabel', canale_eeg, 'XTickLabelRotation', 45);
legend({'Cluster 1','Cluster 2'}, 'Location','northeast');
ylabel('SampEn mediu');
title('SampEn mediu per canal');
grid on;

% Subplot 3 - Silhouette
subplot(1,3,3);
silhouette(X_norm, idx_clustere);
title(sprintf('Silhouette Score = %.3f', mean(s)));

sgtitle('Clustering Nesupervizat K-Means - SampEn EEG GHEATA', ...
        'FontSize',13, 'FontWeight','bold');

cale_folder = 'C:\Users\Bianca\Documents\CercetareIRA\GHEATA\';
saveas(gcf, fullfile(cale_folder, 'CLUSTERING_NESUPERVIZAT_GHEATA.png'));

% 1. Construim tabelul 
TabelClustering = tabel_date(:, {'Fisier'});
TabelClustering.Cluster = idx_clustere;
TabelClustering.Media_Entropie = mean(X, 2);

% 2. cati oameni a pus algoritmul in fiecare grup
nr_oameni_c1 = sum(TabelClustering.Cluster == 1);
nr_oameni_c2 = sum(TabelClustering.Cluster == 2);

TabelClustering.Cluster_Label = repmat({''}, height(TabelClustering), 1);

if nr_oameni_c1 > nr_oameni_c2
    TabelClustering.Cluster_Label(TabelClustering.Cluster == 1) = {'Posibil_Sanatos'};
    TabelClustering.Cluster_Label(TabelClustering.Cluster == 2) = {'Posibil_Bolnav/Atipic'};
else
    TabelClustering.Cluster_Label(TabelClustering.Cluster == 1) = {'Posibil_Bolnav/Atipic'};
    TabelClustering.Cluster_Label(TabelClustering.Cluster == 2) = {'Posibil_Sanatos'};
end

% 4. Salvam fisierul
writetable(TabelClustering, fullfile(cale_folder, 'REZULTATE_CLUSTERING_NESUPERVIZAT4.csv'));
