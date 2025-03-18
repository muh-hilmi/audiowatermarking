function display_results(hasil, logonrz, wtt, M, N, x1)
    hasil.ber0 = mean(abs(wtt - double(logonrz))); % BER sebelum serangan
    hasil.ber1 = mean(abs(wtt - double(logonrz))); % BER setelah serangan (sementara sama karena serangan dimatikan)
    
    logo2d = reshape(logonrz, [M N]);  % Watermark asli
    logo2dt = reshape(wtt, [M N]);     % Watermark yang diekstrak
    
    figure(1), clf
    subplot(121), imshow(uint8(255*logo2d)), title('(a)')
    subplot(122), imshow(uint8(255*logo2dt)), title('(b)')
    
    disp(hasil);
end