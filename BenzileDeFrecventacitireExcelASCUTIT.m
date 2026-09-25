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