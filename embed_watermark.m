function [xw, So] = embed_watermark(x1, logonrz, M, N, B, sub_eksis, mode, alfa, alfass, pn0, pn1, U0, S0, U1, S1, subbad)
    L = B^2;
    xw = zeros(size(x1));
    So(M*N).data = {};
    
    for i = 1:M*N
        % Segmentasi
        seg = alfa * x1(L*i-L+1:L*i);
        
        % Transformasi Wavelet (DWT)
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
            [U, S] = qr(Swdct2);
        end
        
        % Spread Spectrum (mode=1)
        if mode == 1
            So(i).data = {S};
            if logonrz(i) == 1
                Sdgab = S + alfass * S1;
            else
                Sdgab = S + alfass * S0;
            end
            Sisvd = U * Sdgab; % Rekonstruksi QR
        end
        
        % 2D ke 1D
        Sisvd1 = reshape(Sisvd, [(B/2)^2 1]);
        
        % IDST
        if sub_eksis(2) == 1
            Swgab = Sw;
            Swgab(subbad,:) = idst(Sisvd1);
        else
            Swgab(subbad,:) = Sisvd1;
        end
        
        % Rekonstruksi DWT
        cA = idwt(Swgab(1,:), Swgab(2,:), 'haar');
        cD = idwt(Swgab(3,:), Swgab(4,:), 'haar');
        st = idwt(cA, cD, 'haar');
        
        % Gabungkan segmen
        xw(L*i-L+1:L*i) = st / alfa;
    end
end