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