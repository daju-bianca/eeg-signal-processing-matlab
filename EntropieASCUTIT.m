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
writetable(TabelSampEn, fullfile(cale_folder, 'REZULTATE_Entropie_ASCUTIT.csv'));
disp(TabelSampEn);