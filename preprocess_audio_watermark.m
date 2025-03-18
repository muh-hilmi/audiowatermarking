function [x1, fs, logonrz, logobw1d] = preprocess_audio_watermark(x, fs, folderwatermark, M, N, B)
    % Input:
    % x: Sinyal audio dari audiowatermarking1
    % fs: Frekuensi sampling dari audiowatermarking1
    % folderwatermark: Folder tempat watermark berada
    % M, N: Ukuran watermark
    % B: Ukuran segmen
    
    % Resample jika tidak 44100 Hz
    if fs ~= 44100
        x1 = resample(x(:,1), 44100, fs);
        x1 = x1 / (max(abs(x1))*2); % Normalisasi agar maksimum amplitudo jadi 1
        audiowrite("resample.wav", x1, 44100, 'BitsPerSample', 16);
        [x1, fs] = audioread('resample.wav');
    else
        x1 = x(:,1);
    end
    
    % Batasi panjang audio sesuai watermark
    L = B^2;
    x1 = x1(1:L*M*N);
    
    % Membaca dan proses watermark
    logo = imread([folderwatermark 'gambar2.jpeg']);
    logor = imresize(logo(1:round(size(logo,1)*0.75),:,:), [M N]);
    logogray = rgb2gray(logor);
    logobw = ~imbinarize(logogray, 0.0000001);
    logobw1d = double(logobw(:));
    logonrz = 2*logobw1d - 1; % Ubah ke NRZ (-1 dan 1)
end