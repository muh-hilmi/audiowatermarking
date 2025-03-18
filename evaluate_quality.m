function hasil = evaluate_quality(x1, xw, fs, nbit, B)
    hasil.snr = 10*log10((mean(x1.^2))/(mean((x1-xw).^2)));
    hasil.odg = hitungodgsegmen(1, xw, x1, fs, nbit);
    hasil.payload = fs / (B^2);
    
    xw = xw / (max(abs(xw))*2); % Normalisasi agar maksimum amplitudo jadi 1
    audiowrite('temphost.wav', xw, fs, 'BitsPerSample', nbit);
end