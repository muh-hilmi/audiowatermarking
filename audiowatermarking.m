function audiowatermarking1
    % Inisialisasi Parameter Dasar
    M = 64; N = 64;
    Ldek = 3;
    jenis = 2;
    nbit = 16;
    nb = 5;
    mode = 1;
    alfa = 0.1;
    sub_eksis = [2 1 0 2 1]; % DWT, DST, QR, SS
    serangan = [13 64]; % uji coba 0 0 dan 13 64
    
    folderhost = [pwd '/host_audio/'];
    folderwatermark = [pwd '/watermark/'];
    
    % Audio input
    [x, fs] = audioread([folderhost 'africa-toto.wav']);
    
    % Parameter uji coba
    subband_range = [1, 2, 3, 4]; % Pilihan subband
    B_range = [2, 4, 6, 8, 16, 32]; % Ukuran segmen
    alfass_range = 0.001:0.05:0.1; % Kekuatan SS
    
    % Inisialisasi hasil
    results = struct('subband', [], 'alfass', [], 'B', [], 'snr', [], 'odg', [], 'payload', [], 'ber0', [], 'ber1', []);
    test_idx = 1;
    
    % Loop subband
    for subband = subband_range
        % Loop B
        for B = B_range
            stop_alfass = false; % Flag untuk menghentikan loop alfass jika SNR < 15
            
            % Loop alfass
            for i_alfass = 1:length(alfass_range)
                if stop_alfass
                    break; % Jika sudah berhenti, lanjut ke B berikutnya
                end
                
                alfass = alfass_range(i_alfass);
                
                try
                    disp(['Uji coba ', num2str(test_idx), ...
                          ': subband=', num2str(subband), ', B=', num2str(B), ', alfass=', num2str(alfass)]);
                    
                    % Preprocessing
                    [x1, fs, logonrz, logobw1d] = preprocess_audio_watermark(x, fs, folderwatermark, M, N, B);
                    % Inisialisasi SS
                    [pn0, pn1, U0, S0, U1, S1] = init_spread_spectrum(sub_eksis, B);
                    % Embedding
                    [xw, So] = embed_watermark(x1, logonrz, M, N, B, sub_eksis, mode, alfa, alfass, pn0, pn1, U0, S0, U1, S1, subband);
                    % Pengujian
                    hasil = evaluate_quality(x1, xw, fs, nbit, B);
                    
                    % Jika SNR di bawah 15, hentikan loop alfass
                    if hasil.snr < 15
                        disp(['SNR = ', num2str(hasil.snr), ' dB (di bawah 15), hentikan iterasi alfass...']);
                        stop_alfass = true;
                        break;
                    end
                    
                    % Serangan
                    xwn = apply_attack(xw, fs, serangan, nbit);
                    % Ekstraksi
                    wtt = extract_watermark(xwn, M, N, B, sub_eksis, mode, alfa, alfass, So, U0, S0, U1, S1, pn0, pn1, subband);
                    
                    % Hitung BER
                    % hasil.ber0 = mean(abs(wt - double(logobw1d))); % BER sebelum serangan
                    hasil.ber1 = mean(abs(wtt - double(logobw1d))); % BER setelah serangan
                    
                
                    % Simpan hasil
                    results(test_idx).subband = subband;
                    results(test_idx).alfass = alfass;
                    results(test_idx).B = B;
                    results(test_idx).snr = hasil.snr;
                    results(test_idx).odg = hasil.odg;
                    results(test_idx).payload = hasil.payload;
                    results(test_idx).ber0 = hasil.ber0;
                    results(test_idx).ber1 = hasil.ber1;
                    
                    test_idx = test_idx + 1;
                catch ME
                    disp(['Error pada subband=', num2str(subband), ', B=', num2str(B), ', alfass=', num2str(alfass)]);
                    disp(ME.message);
                    % Simpan hasil kosong untuk uji coba yang gagal
                    results(test_idx).subband = subband;
                    results(test_idx).alfass = alfass;
                    results(test_idx).B = B;
                    results(test_idx).snr = NaN;
                    results(test_idx).odg = NaN;
                    results(test_idx).payload = NaN;
                    results(test_idx).ber0 = NaN;
                    results(test_idx).ber1 = NaN;
                    test_idx = test_idx + 1;
                end
            end
        end
    end
    
    % Tampilkan semua hasil
    disp('Hasil Uji Coba:');
    disp(struct2table(results));
    
    % Simpan ke file MAT
    save('watermarking_results.mat', 'results');
    
    % Simpan ke file CSV
    results_table = struct2table(results);
    writetable(results_table, 'watermarking_results.csv');
    disp('Hasil telah disimpan ke watermarking_results.csv');
end
