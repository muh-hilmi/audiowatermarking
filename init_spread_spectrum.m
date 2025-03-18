function [pn0, pn1, U0, S0, U1, S1] = init_spread_spectrum(sub_eksis, B)
    % Generate pncode berdasarkan wavelet
    if sub_eksis(1) == 1 || sub_eksis(1) == 0 % SWT atau tanpa wavelet
        pn0 = 2*randi([0 1], B, B) - 1;
        pn1 = 2*randi([0 1], B, B) - 1;
    elseif sub_eksis(1) == 2 % DWT
        pn0 = 2*randi([0 1], B/2, B/2) - 1;
        pn1 = 2*randi([0 1], B/2, B/2) - 1;
    end
    
    % Dekomposisi untuk SS
    if sub_eksis(4) == 1 % SVD
        [U0, S0, V0] = svd(pn0);
        [U1, S1, V1] = svd(pn1);
    elseif sub_eksis(4) == 2 % QR
        [U0, S0] = qr(pn0);
        [U1, S1] = qr(pn1);
    end
end