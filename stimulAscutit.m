% procesarea unui singur set de date 
data = readtable('C:\Users\Bianca\Documents\CercetareIRA\Stimul_ascutit (1)\BBT-E12-AAB038-2024-11-11_09-52-35\EEG.csv');
eeg_data = table2array(data(:, 5:16))';
% verific dimensiunea
disp(size(eeg_data))

%% Procesare perosane-STIMUL ASCUTIT
cd_proiect = 'C:\Users\Bianca\Documents\CercetareIRA\Stimul_ascutit (1)';
datasetList = dir(fullfile(cd_proiect, 'BBT*')); 
datasetList = datasetList([datasetList.isdir]); 

fprintf('Sau găsit %d persoane...\n', length(datasetList));

for i = 1:length(datasetList)
    dataName = datasetList(i).name;
    fprintf('\n>>> Procesăm persoana %d/%d: %s <<<\n', i, length(datasetList), dataName);
    
    % 2. Importul datelor (.CSV)
    cale_csv = fullfile(cd_proiect, dataName, 'EEG.csv');
    temp_table = readtable(cale_csv);
    eeg_data = temp_table{:, 5:16}'; 
    
    % Import în EEGLAB
    EEG = pop_importdata('dataformat','array','nbchan',12,'data','eeg_data','srate',256);
    EEG.setname = dataName;

    % 3. FILTRU TRECE-SUS (1 Hz)
    EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'plotfreqz', 0);

    % 4. FILTRU TRECE-JOS (40 Hz)  
    EEG = pop_eegfiltnew(EEG, 'hicutoff', 40, 'plotfreqz', 0);

    % 5. FILTRU NOTCH (50 Hz) 
    EEG = pop_eegfiltnew(EEG, 'locutoff', 48, 'hicutoff', 52, 'revfilt', 1, 'plotfreqz', 0);

    % 6. CURĂȚARE SPIKE 
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 5, 'ChannelCriterion', 0.8, ...
        'BurstCriterion', 20, 'WindowCriterion', 'off', 'BurstRejection', 'off');

    % 7. RE-REFERENȚIERE 
    EEG = pop_reref(EEG, []);

    % 8. SALVARE AUTOMATĂ (.set)
    nume_salvare = [dataName '_Procesat.set'];
    EEG = pop_saveset(EEG, 'filename', nume_salvare, 'filepath', cd_proiect);
    ALLEEG = []; EEG = [];
end
%% Analiza benzilor GHEATA - Stimul Ascutit
cale_folder = 'C:\Users\Bianca\Documents\CercetareIRA\GHEATA\'; 
lista_fisiere = dir(fullfile(cale_folder, '*Procesat.set'));
nume_electrozi = {'Fp1', 'Fp2', 'F7', 'F3', 'F4', 'F8', 'P7', 'P3', 'P4', 'P8', 'O1', 'O2'};

% Inițializare EEGLAB pentru a curăța memoria
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
rezultate_matrice = zeros(length(lista_fisiere), 5);
nume_subiecti = cell(length(lista_fisiere), 1);
fprintf('Procesare GHEATA pentru %d subiecți...\n', length(lista_fisiere));

for i = 1:length(lista_fisiere)
    EEG = pop_loadset('filename', lista_fisiere(i).name, 'filepath', cale_folder);
    for ch = 1:EEG.nbchan
        EEG.chanlocs(ch).labels = nume_electrozi{ch};
    end
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

% 6. Salvare Rezultate
TabelGheata = array2table(rezultate_matrice, 'VariableNames', {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'});
TabelGheata.Fisier = nume_subiecti;
TabelGheata = TabelGheata(:, [6 1 2 3 4 5]); 
writetable(TabelGheata, fullfile(cale_folder, 'REZULTATE_FINALE_GHEATA.csv'));
eeglab redraw; 
%%  Afisarea benzilor de frecventa
clear; clc; close all;
[fisier, cale] = uigetfile('*.set', 'Selecteaza fisierul EEG (Stimul Ascutit)');
if isequal(fisier, 0)
    disp('eroare');
    return;
end

EEG = pop_loadset('filename', fisier, 'filepath', cale);
data = EEG.data(1, :); % Extragem datele de pe primul canal (Fp1 de obicei)

Fs = EEG.srate;        % Frecvența de eșantionare preluată din aparat
L = length(data);      % Lungimea semnalului
t = (0:L-1)*(1/Fs);    % Vectorul de timp
timp_afisare = 5;      % Fereastra de vizualizare (secunde)
NFFT = 2^nextpow2(L);  % Optimizare pentru FFT
f = Fs/2*linspace(0, 1, NFFT/2); % Vectorul de frecvențe

% --- BANDA DELTA (0.5 - 4 Hz) ---
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
title('Forma de undă pentru banda DELTA (ASCUTIT)');
xlabel('Timp (secunde)'); ylabel('Amplitudine (\muV)');
xlim([0 timp_afisare]); 

Y1 = fft(x1, NFFT)/L;
h8 = figure; 
plot(f, 2*abs(Y1(1:NFFT/2)), 'r'); 
title('Spectrul de amplitudine unilateral pentru banda DELTA (ASCUTIT)');
xlabel('Frecvență (Hz)'); ylabel('|Amplitudine|');
xlim([0 40]); 

% --- BANDA THETA (4-7 Hz) ---
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
title('Forma de undă pentru banda THETA (ASCUTIT)');
xlabel('Timp (secunde)'); ylabel('Amplitudine (\muV)');
xlim([0 timp_afisare]); 

Y2 = fft(x2, NFFT)/L;
h9 = figure; 
plot(f, 2*abs(Y2(1:NFFT/2)), 'm'); 
title('Spectrul de amplitudine unilateral pentru banda THETA (ASCUTIT)');
xlabel('Frecvență (Hz)'); ylabel('|Amplitudine|');
xlim([0 40]); 

% --- BANDA ALPHA (8-12 Hz) ---
Fstop1 = 7.5;           
Fpass1 = 8;             
Fpass2 = 12;            
Fstop2 = 12.5;          

[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
b3 = firpm(N, Fo, Ao, W, {dens});
Hd3 = dfilt.dffir(b3);
x3 = filter(Hd3, data);

h6 = figure; plot(t, x3, 'g'); 
title('Forma de undă pentru banda ALPHA (ASCUTIT)');
xlabel('Timp (secunde)'); ylabel('Amplitudine (\muV)');
xlim([0 timp_afisare]); 

Y3 = fft(x3, NFFT)/L;
h10 = figure; 
plot(f, 2*abs(Y3(1:NFFT/2)), 'g'); 
title('Spectrul de amplitudine unilateral pentru banda ALPHA (ASCUTIT)');
xlabel('Frecvență (Hz)'); ylabel('|Amplitudine|');
xlim([0 40]); 

% --- BANDA BETA (12-30 Hz) ---
Fstop1 = 11.5;          
Fpass1 = 12;            
Fpass2 = 30;            
Fstop2 = 30.5;          

[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
b4 = firpm(N, Fo, Ao, W, {dens});
Hd4 = dfilt.dffir(b4);
x4 = filter(Hd4, data);

h7 = figure; plot(t, x4, 'b'); 
title('Forma de undă pentru banda BETA (ASCUTIT)');
xlabel('Timp (secunde)'); ylabel('Amplitudine (\muV)');
xlim([0 timp_afisare]); 

Y4 = fft(x4, NFFT)/L;
h11 = figure; 
plot(f, 2*abs(Y4(1:NFFT/2)), 'b'); 
title('Spectrul de amplitudine unilateral pentru banda BETA (ASCUTIT)');
xlabel('Frecvență (Hz)'); ylabel('|Amplitudine|');
xlim([0 40]);
%% Cod entropie 
%  PhysioNet: https://physionet.org/content/sampen/1.0.0/
cale_folder = 'C:\Users\Bianca\Documents\CercetareIRA\Stimul_ascutit (1)';
lista_fisiere = dir(fullfile(cale_folder, '*.set'));
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
    for ch = 1:EEG.nbchan
        semnal = double(EEG.data(ch, 1:min(3000, EEG.pnts)));
        e = sampen(semnal, m, r, 'chebychev');
        rezultate_sampen(i, ch) = e;
        
    end

    disp(['Subiect ', num2str(i), ' din ', num2str(length(lista_fisiere)), ' finalizat.']);
end

TabelSampEn = array2table(rezultate_sampen, 'VariableNames', nume_electrozi);
TabelSampEn.Fisier = nume_subiecti;
TabelSampEn = TabelSampEn(:, [13 1:12]);
writetable(TabelSampEn, fullfile(cale_folder, 'REZULTATE_SAMPEN_ASCUTIT.csv'));
disp(TabelSampEn);
%% CLUSTERING NESUPERVIZAT - K-MEANS EEG ASCUTIT
% Sursa 1: MathWorks. Cluster Analysis Example
%Sursa 2:AI
cale_folder_ascutit = 'C:\Users\Bianca\Documents\CercetareIRA\Stimul_ascutit (1)\';
fisier_date = fullfile(cale_folder_ascutit, 'REZULTATE_SAMPEN_ASCUTIT.csv');
tabel_date = readtable(fisier_date);

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
figure('Name','Elbow Method - ASCUTIT','Position',[100 100 500 400]);
sumd_vec = zeros(1,6);
for k = 1:6
    [~,~,sumd] = kmeans(X_norm, k, 'Replicates', 10, 'Display', 'off');
    sumd_vec(k) = sum(sumd);
end
plot(1:6, sumd_vec, '-o', 'LineWidth', 2, 'MarkerFaceColor','b', 'MarkerSize', 8);
xlabel('Numar clustere (k)');
ylabel('Suma distante intra-cluster');
title('Elbow Method - Numar optim clustere (ASCUTIT)');
xline(2,'--r','k=2','LineWidth',1.5);
grid on;

% 5. K-MEANS k=2
rng('default');  
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

% 7. PURITATE CLUSTERE
media_c1 = mean(X(idx_clustere==1, :));
media_c2 = mean(X(idx_clustere==2, :));
fprintf('\nMedia SampEn per cluster:\n');
fprintf('%-6s | %-12s | %-12s | %s\n', 'Canal','Cluster 1','Cluster 2','Diferenta');
fprintf('%s\n', repmat('-',1,45));
for j = 1:length(canale_eeg)
    fprintf('%-6s | %-12.4f | %-12.4f | %.4f\n', ...
            canale_eeg{j}, media_c1(j), media_c2(j), media_c1(j)-media_c2(j));
end
figure('Name','Clustering Nesupervizat - ASCUTIT', ...
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
title('Clustere K-Means pe PCA (ASCUTIT)');
legend({'Cluster 1','Cluster 2','Centroizi'}, 'Location','best','FontSize',8);
grid on;

% Subplot 2 - Media SampEn 
subplot(1,3,2);
bar([media_c1; media_c2]', 'grouped');
set(gca,'XTickLabel', canale_eeg, 'XTickLabelRotation', 45);
legend({'Cluster 1','Cluster 2'}, 'Location','northeast');
ylabel('SampEn mediu');
title('SampEn mediu per canal (ASCUTIT)');
grid on;

% Subplot 3 - Silhouette
subplot(1,3,3);
silhouette(X_norm, idx_clustere);
title(sprintf('Silhouette Score = %.3f', mean(s)));
sgtitle('Clustering Nesupervizat K-Means - SampEn EEG ASCUTIT', ...
        'FontSize',13, 'FontWeight','bold');
saveas(gcf, fullfile(cale_folder_ascutit, 'CLUSTERING_NESUPERVIZAT_ASCUTIT.png'));
SampEn_Global = mean(X, 2); 

TabelClustering = tabel_date(:, {'Fisier'});
TabelClustering.SampEn_Mediu = SampEn_Global; % Adăugăm media SampEn
TabelClustering.Cluster = idx_clustere;
TabelClustering.Cluster_Label = repmat({''}, height(TabelClustering), 1);
TabelClustering.Cluster_Label(idx_clustere==1) = {'Posibil_Sanatos'};
TabelClustering.Cluster_Label(idx_clustere==2) = {'Posibil_Bolnav'};
% Salvăm tabelul în format CSV
writetable(TabelClustering, fullfile(cale_folder_ascutit, 'REZULTATE_CLUSTERING_NESUPERVIZAT_ASCUTIT1.csv'));
fprintf('\n=== GATA! Salvat in folderul Stimul_ascutit (1) ===\n');
disp(TabelClustering);