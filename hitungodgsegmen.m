function odg=hitungodgsegmen(odg_cal,wavt,wav,fs,bit)
tipeodg='baru';

    audiowrite('aw_temp.wav',wavt,fs,'BitsPerSample',bit)
    audiowrite('aw_host.wav',wav,fs,'BitsPerSample',bit)


if odg_cal==1
%     if length(wav)<1200
%         NN=1200-length(wav);
%         wav=[wav ones(1,NN)];
%         wavt=[wavt ones(1,NN)];
%         audiowrite('aw_temp.wav',wavt,fs,'BitsPerSample',bit)
%         audiowrite('aw_host.wav',wav,fs,'BitsPerSample',bit)
%         odg=PQevalAudio ('aw_host.wav','aw_temp.wav');
%     else
        if strcmp(tipeodg,'lama')
            odg=PQevalAudio ('aw_host.wav','aw_temp.wav');
        elseif strcmp(tipeodg,'baru')
%             [x teksodg]=system(['eaqual -srate 44100 -fref "' file_host_cut '" -ftest "' file_aw '"']);
            [x teksodg]=system(['eaqual -srate 44100 -fref "aw_host.wav" -ftest "aw_temp.wav"']);
            posodga=findstr(teksodg,'ODG');
            pos9=find(teksodg(posodga:posodga+7)==9);
            pos10=find(teksodg(posodga:posodga+20)==10);
            odg=str2num(teksodg(posodga+pos9:posodga+pos10-1));
        end
%     end
else
    odg=-4;
end