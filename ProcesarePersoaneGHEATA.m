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