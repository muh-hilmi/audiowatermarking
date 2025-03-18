function wtt = extract_watermark(xwn, M, N, B, sub_eksis, mode, alfa, alfass, So, U0, S0, U1, S1, pn0, pn1, subbad)
    L = B^2;
    wtt = zeros(M*N, 1);
    
    for i = 1:M*N
        % Segmentasi
        seg = alfa * xwn(L*i-L+1:L*i);
        
        % DWT
        if sub_eksis(1) == 2
            [A, D] = dwt(seg', 'haar');
            [LL, LH] = dwt(A, 'haar');
            [HL, HH] = dwt(D, 'haar');
            Sw = [LL; LH; HL; HH];
        end
        
        % DST
        if sub_eksis(2) == 1
            sdst = dst(Sw(subbad,:));
        else
            sdst = Sw(subbad,:);
        end
        
        % 1D ke 2D
        Swdct2 = reshape(sdst, [B/2 B/2]);
        
        % QR Decomposition
        if sub_eksis(4) == 2
            [~, S] = qr(Swdct2);
        end
        
        % Ekstraksi SS
        if mode == 1
            X0 = mean(mean(abs((U0*(S-cell2mat(So(i).data))/alfass).*pn0))); % pn0 digunakan di sini
            X1 = mean(mean(abs((U1*(S-cell2mat(So(i).data))/alfass).*pn1))); % pn1 digunakan di sini
            if X0 >= X1
                wtt(i,1) = 0;
            else
                wtt(i,1) = 1;
            end
        end
    end
end