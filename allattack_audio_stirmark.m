function allattack_audio_stirmark(serangan,Filewater,Fileattack,nbit)
matlabver='baru';
% folder=[folder '\'];
% folderserang=[folderserang '\'];
% if isempty(dir(folderserang))
%     mkdir(folderserang)
% end
% folderserang
jenis=serangan(1,1);
indserangan=serangan(1,2);
if length(jenis)==3
    indserangan=serangan(1,2:end);
end
A=1;
folderattack=[pwd '\'];
% dire=dir(folder);
% h=waitbar(0,'Wait ....');
dump=1;
% for i=1:size(dire,1)
for i=1:1
    %     waterfile=[folder dire(i).name];
    waterfile=Filewater;
    %     attackfile=[folderserang num2str(serangan(1,1)) num2str(serangan(1,2))  '-' dire(i).name(3:end) '.wav'];
    attackfile=Fileattack;
    % p=find(attackfile=='\');
    % folderattack=Fileattack(1:p(end));
    %     if ~dire(i).isdir
    %         switch dire(i).name(1,end-2:end)
    switch attackfile(1,end-2:end)
        case 'wav'
                [wave,Fs]=audioread(waterfile);
            wave=wave(:,1);
            if jenis==0 %Tanpa serangan
                    audiowrite(attackfile,A*wave,Fs,'BitsPerSample',nbit)
            elseif jenis==1%LPF fleksibel
                f=indserangan;
                %                     f=[ 1000];
                Fs=44100;
                [b,a]=butter(2,f/(Fs/2));
                wave2=filtfilt(b,a,wave);
                %                     wavwrite(A*wave2,Fs,nbit,attackfile)
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==2 %BPF
                if indserangan==1
                    f=[ 100 3000];
                elseif indserangan==2
                    f=[ 100 6000];
                elseif indserangan==3
                    f=[ 100 9000];
                elseif indserangan==4
                    f=[ 50 6000];
                elseif indserangan==5
                    f=[ 25 6000];
                end
                
                Fs=44100;
                [b,a]=butter(2,f/(Fs/2));
                wave2=filtfilt(b,a,wave);
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==3 %Rekuantisasi
                    audiowrite(attackfile,A*wave,Fs,'BitsPerSample',indserangan)
            elseif jenis==4 %Stereo 2 Mono
                wave2=mean(wave,2);
                %                     wavwrite(A*wave2,Fs,nbit,attackfile)
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==5 % Noise
                %                 if indserangan==1
                %                     SNRdb=0;
                %                 elseif indserangan==2
                %                     SNRdb=10;
                %                 elseif indserangan==3
                %                     SNRdb=20;
                %                 end
                SNRdb=indserangan;
                G=sqrt(10^(-SNRdb/10));
                wavemax=max(abs(wave));
                wavepow=mean(wave.^2,1);
                pinkn=G*sqrt(wavepow)*pinknoise(length(wave)).';
                pinknmax=max(abs(pinkn));
                pinknpow=mean(pinkn.^2,1);
                whiten=G*sqrt(wavepow)*wgn(1,length(wave),0,'real').';
                whitenmax=max(abs(whiten));
                whitenpow=mean(whiten.^2,1);
%                 wave2=wave+pinkn+whiten;
                wave2=wave+whiten;
                %                     wavwrite(A*wave2,Fs,nbit,attackfile)
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==6 %Resampling
                if indserangan==1
                    Fsb=11025;
                elseif indserangan==2
                    Fsb=16000;
                elseif indserangan==3
                    Fsb=22050;
                elseif indserangan==4
                    Fsb=24000;
                elseif indserangan==5
                    Fsb=48000;
                end
                wave2=resample(resample(wave(:,1),Fsb,Fs),Fs,Fsb);
                if length(wave2)>length(wave)
                    wave2=wave2(1:length(wave));
                elseif length(wave2)<length(wave)
                    wave2=[wave2;zeros(length(wave)-length(wave2),1)];
                end
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==7 %Time Scale Modification
                if indserangan==1
                    persen=1;
                elseif indserangan==2
                    persen=2;
                elseif indserangan==3
                    persen=3;
                elseif indserangan==4
                    persen=4;
                end
                var=1+persen/100;vari=1/var;
                %                     var=1/(1-persen/100);vari=1/var;
                %                 system(['sox -q -V1 "',waterfile,'" "',[folderattack 'temp.wav'],'" tempo ' num2str(var)]);
                %                 system(['sox -q -V1 temp.wav temp1.wav tempo ' num2str(vari)]);
                system(['sox -q -V1 "',waterfile,'" "',[folderattack 'temp.wav'],'" tempo ' num2str(var)]);
                system(['sox -q -V1 "',[folderattack 'temp.wav'],'" "',[folderattack 'temp1.wav'],'" tempo ' num2str(vari)]);
                
                if strcmp(matlabver,'lama')
                    [wave2,fs]=wavread([folderattack 'temp1.wav']);
                else
                    [wave2,fs]=audioread([folderattack 'temp1.wav']);
                end
                if length(wave2)>length(wave)
                    wave2=wave2(1:length(wave));
                elseif length(wave2)<length(wave)
                    wave2=[wave2;zeros(length(wave)-length(wave2),1)];
                end
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==8 %Speed Change
                if indserangan==1
                    persen=1;
                elseif indserangan==2
                    persen=5;
                elseif indserangan==3
                    persen=10;
                end
                var=1+persen/100;vari=1/var;
                system(['sox -q -V1 "',waterfile,'" "',[folderattack 'temp.wav'],'" speed ' num2str(var)]);
                system(['sox -q -V1 "',[folderattack 'temp.wav'],'" "',[folderattack 'temp1.wav'],'" speed ' num2str(vari)]);
                
                    [wave2,fs]=audioread([folderattack 'temp1.wav']);
                
                if length(wave2)>length(wave)
                    wave2=wave2(1:length(wave));
                elseif length(wave2)<length(wave)
                    wave2=[wave2;zeros(length(wave)-length(wave2),1)];
                end
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
                
                
            elseif jenis==9 %Pitch Shifting
                if indserangan==1
                    persen=1;
                elseif indserangan==2
                    persen=2;
                elseif indserangan==3
                    persen=3;
                elseif indserangan==4
                    persen=4;
                end
                var=persen/100*1200;
                system(['sox -q -V1 "',waterfile,'" "',attackfile,'" pitch ' num2str(var)]);
            elseif jenis==10 %
                wave2=equalizer(wave);
                %                     wavwrite(A*wave2,Fs,nbit,attackfile)
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==11 %Gema/Echo
                system(['sox -q -V1 "',waterfile,'" "',[folderattack 'temp1.wav'],'" echo 1 0.3 100 ' num2str(A)]);
                
                    [wave2,fs]=audioread([folderattack 'temp1.wav']);
                if length(wave2)>length(wave)
                    wave2=wave2(1:length(wave));
                elseif length(wave2)<length(wave)
                    wave2=[wave2;zeros(length(wave)-length(wave2),1)];
                end
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==12 %Cropping
                if indserangan==1
                    persen=25;
                elseif indserangan==2
                    persen=50;
                elseif indserangan==3
                    persen=75;
                elseif indserangan==4
                    persen=100;
                end
                wave2=[wave(floor(persen/100*Fs)+1:end,1) ;zeros(floor(persen/100*Fs),1)];
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==13 %MP3
                rate=indserangan;
                
%                 system(['sox -q -V1 -C 32 "', [waterfile], '" "', ['temp1.mp3'], '" -C ' num2str(rate)]);
%    system(['sox -q -V1 temphost.wav -C 32
%        temp1.mp3']);
%                 
%                 system(['sox -q -V1 "', [folderattack waterfile], '" "', [folderattack 'temp1.mp3'], '" -C ' num2str(rate)]);
%    
%                 system(['sox -q -V1 "', [folderattack 'temphost.wav'], ...
%         '" "', [folderattack 'temp1.mp3'], '" -C ' num2str(rate)]);
%     
%                 system(['sox temp1.mp3 temp1.wav']);
% 
%                 
%                 inputFile = [folderattack 'temphost.wav'];
% outputFile = [folderattack 'temp1.mp3'];
% rate = 32;  % Tingkat kompresi MP3 dalam kbps
% 
% % Perintah SoX
% cmd = ['sox -q -V1 "' inputFile '" "' outputFile '" -C ' num2str(rate)];
% 
% % Menjalankan perintah SoX
% status = system(cmd);
                
                
%                 system(['sox -q -V1 "',waterfile,'" -C "' num2str(rate) ,[folderattack 'temp.mp3'],'" speed ' num2str(rate)]);
%                 system(['sox -q -V1 "',[folderattack 'temp.wav'],'" "',[folderattack 'temp1.wav'],'" speed ' num2str(vari)]);
                
%                     [wave2,fs]=audioread([folderattack 'temp1.wav']);

                
                system(['ffmpeg -y -i "',waterfile,'" -loglevel error -b:a ' num2str(rate) 'k "',[folderattack 'temp.mp3'],'"']);
                system(['ffmpeg -y -i "',[folderattack 'temp.mp3'],'" -loglevel error "',attackfile '"']);
                
                    [wave2,fs]=audioread(attackfile);
                if length(wave2)>length(wave)
                    wave2=wave2(1:length(wave));
                elseif length(wave2)<length(wave)
                    wave2=[wave2;zeros(length(wave)-length(wave2),1)];
                end
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==14 %AAC
                if indserangan==1
                    rate=32;
                elseif indserangan==2
                    rate=64;
                elseif indserangan==3
                    rate=128;
                elseif indserangan==4
                    rate=192;
                end
                [p,q]=system(['ffmpeg -y -i "',waterfile,'" -loglevel error -strict experimental  -b ' num2str(rate) 'k "',[folderattack 'temp.aac'],'"']);
                [p,q]=system(['ffmpeg -y -i "',[folderattack 'temp.aac'],'" -loglevel error "',attackfile '"']);
                    [wave2,fs]=audioread(attackfile);
                if length(wave2)>length(wave)
                    wave2=wave2(1:length(wave));
                elseif length(wave2)<length(wave)
                    wave2=[wave2;zeros(length(wave)-length(wave2),1)];
                end
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==15 %MP4
                if indserangan==1
                    rate=32;
                elseif indserangan==2
                    rate=64;
                elseif indserangan==3
                    rate=128;
                elseif indserangan==4
                    rate=192;
                end
                [p,q]=system(['ffmpeg -y -i "',waterfile,'" -loglevel error  -strict -2 -b ' num2str(rate) 'k "',[folderattack 'temp.mp4'],'"']);
                [p,q]=system(['ffmpeg -y -i "',[folderattack 'temp.mp4'],'" -loglevel error "',attackfile,'"']);
                    [wave2,fs]=audioread(attackfile);
                if length(wave2)>length(wave)
                    wave2=wave2(1:length(wave));
                elseif length(wave2)<length(wave)
                    wave2=[wave2;zeros(length(wave)-length(wave2),1)];
                end
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==16 %Dolby
                if indserangan==1
                    rate=32;
                elseif indserangan==2
                    rate=64;
                elseif indserangan==3
                    rate=128;
                elseif indserangan==4
                    rate=192;
                end
                [p,q]=system(['ffmpeg -y -i "',waterfile,'" -loglevel error  -b ' num2str(rate) 'k "',[folderattack 'temp.ac3'],'"']);
                [p,q]=system(['ffmpeg -y -i "',[folderattack 'temp.ac3'],'" -loglevel error "',attackfile '"']);
                    [wave2,fs]=audioread(attackfile);
                if length(wave2)>length(wave)
                    wave2=wave2(1:length(wave));
                elseif length(wave2)<length(wave)
                    wave2=[wave2;zeros(length(wave)-length(wave2),1)];
                end
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==17 %Delay
                wave2=[zeros(indserangan,1);wave];
                    audiowrite(attackfile,A*wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==18 %AddBrumm
                t=linspace(0,1,length(wave)).';
                A=indserangan;
                f=55;
                wave2=wave+A*sin(2*pi*f*t);
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==19 %AddSinus
                t=linspace(0,1,length(wave)).';
                A=indserangan;
                f=3000;
                %                 f=indserangan(1,2);
                wave2=wave+A*sin(2*pi*f*t);
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==20 %AddNoise
                
                PN=randn(size(wave));
                As=indserangan;
                b=1;
                k=mod(PN,As);
                a=1-As/A;
                wave2=a*wave+b*k;
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==21 %AddDynNoise
                PN=randn(size(wave));
                As=indserangan;
                b=1;a=1;
                k=mod(PN,abs(wave.*As/100)+1);
                wave2=a*wave+b*k;
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==22 %AddFFTNoise
                
                PN=randn(size(wave));
                PN=PN(1:end/2).';
                As=indserangan;
                b=1;
                k=mod(PN,As);
                a=1-As/A;
                Wave=fft(wave).';
                Wave=Wave(:,1:end/2);
                
                wave2=a*Wave+b*k;
                
                wave2(:,end+2:length(wave2)*2) = real(wave2(:,end:-1:2)) -1i*imag(wave2(:,end:-1:2));
                wavet=real(ifft(wave2.').');
                
                
                    audiowrite(attackfile,wavet,Fs,'BitsPerSample',nbit)
            elseif jenis==23 %Denoise
                [wave2,cfs] = cmddenoise(wave,'sym4',5,'s');
                
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==24 %LSBzero
                wave1=int16((2^15-1)*wave);
                wave3= wave1- rem(wave1, 2);
                wave2=double(wave3)/(2^15-1);
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==25 %Echo
                delay=indserangan;
                a=2;b=2;
                wave2=a*wave+b*[zeros(delay,1);wave(1:end-delay,1)];
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==26 %Amplify
                As=indserangan;
                wave2=wave*As/100;
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==27 %Normalizer1 time domain
                wave2=indserangan/2^15*wave/max(abs(wave));
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==28 %Normalizer2 freq. domain
                wavef=fft(wave);
                wavefn=indserangan/2^15*wavef/max(abs(wavef));
                wave2=ifft(wavefn);
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==29 %Compressor
                CV=2.1;thrdb=-6.1;
                wave2=(10*log10(wave.^2)>thrdb)*CV+(10*log10(wave.^2)<=thrdb).*wave;
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==30 %Bass Boost
                gaindb=indserangan;
                ft=2^7;
                orde=3;
                gain=10.^(gaindb/20);
                fs=44100;
                fc=ft/(fs/2);
                [num1,den1]=butter(orde,fc);
                wave2=wave+gain*filtfilt(num1,den1,wave);
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==31 %RC High-Pass
                thr=150;fs=44100;
                B=exp(-2*pi*thr/fs);
                A0=(1+B)/2;A1=-(1+B)/2;
                
                wave2=wave;
                wave2(1)=0;
                for ii=2:length(wave)
                    wave2(ii)=A0*wave(ii)+A1*wave(ii-1)+B*wave2(ii-1);
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==32 %RC Low-Pass
                thr=15000;fs=44100;
                B=exp(-2*pi*thr/fs);
                A=(1-B);
                
                wave2=wave;
                wave2(1)=0;
                for ii=2:length(wave)
                    wave2(ii)=A*wave(ii)+B*wave2(ii-1);
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==33 %FFT_HLPassQuick
                wavef=fft(wave);
                
                f1=100;
                f2=16000;
                Wavef=wavef(1:end/2).';
                F1=round(f1/(Fs/2)*length(Wavef));
                F2=round(f2/(Fs/2)*length(Wavef));
                Wavez=Wavef;
                Wavez([1:F1 F2:end])=0;
                Wavez(:,end+2:length(Wavez)*2) = real(Wavez(:,end:-1:2)) -1i*imag(Wavez(:,end:-1:2));
                wave2=real(ifft(Wavez.').');
                
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==34 %Stat1
                
                wave2=wave;
                for ii=1:length(wave)-1
                    wave2(ii)=(wave(ii)+wave(ii+1))/2;
                end
                wave2(end)=wave(end);
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==35 %Stat2
                
                wave2=wave;
                for ii=2:length(wave)-1
                    wave2(ii)=(wave(ii-1)+3*wave(ii)+wave(ii+1))/5;
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==36 %FFTStat1
                wavef=fft(wave);
                Wavef=wavef(1:end/2).';
                
                wave1=Wavef;
                for ii=1:length(Wavef)-1
                    wave1(ii)=(Wavef(ii)+Wavef(ii+1))/2;
                end
                %                 wave1(end)=Wavef(end);
                
                
                wave1(:,end+2:length(wave1)*2) = real(wave1(:,end:-1:2)) -1i*imag(wave1(:,end:-1:2));
                wave2=real(ifft(wave1.').');
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==37 %Smooth1
                
                wave2=wave;
                for ii=2:length(wave)-1
                    if wave(ii-1)>wave(ii) && wave(ii+1)>wave(ii)
                        wave2(ii)=(wave(ii-1)+wave(ii+1))/2;
                    elseif wave(ii-1)<wave(ii) && wave(ii+1)<wave(ii)
                        wave2(ii)=(wave(ii-1)+wave(ii+1))/2;
                    end
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==38 %Smooth2
                
                wave2=wave;
                for ii=2:length(wave)-1
                    if wave(ii-1)<wave(ii) && wave(ii+1)>wave(ii)
                        wave2(ii)=(wave(ii-1)+wave(ii+1))/2;
                    elseif wave(ii-1)>wave(ii) && wave(ii+1)<wave(ii)
                        wave2(ii)=(wave(ii-1)+wave(ii+1))/2;
                    end
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==39 %Invert
                
                wave2=-wave;
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==40 %FFTInvert
                
                wavef=fft(wave);
                wave1=-wavef(1:end/2).';
                
                wave1(:,end+2:length(wave1)*2) = real(wave1(:,end:-1:2)) -1i*imag(wave1(:,end:-1:2));
                wave2=real(ifft(wave1.').');
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
            elseif jenis==41 %CopySample belum selesai
                DI=6000;C=2000;Period=10000;
                
                wave2=wave;
                for ii=1:length(wave)
                    ii
                    z=mod(ii,Period+C);
                    z
                    if z>=0 && z<DI
                        wave2(ii)=wave(z);
                    elseif z>=DI && z<DI+C
                        wave2(ii)=wave(z-DI);
                    elseif z>=DI+C && z<Period+C
                        wave2(ii)=wave(z-C);
                    end
                end
                %                 wave2(1)=wave(1);
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==42 %FlippSample
                DI=6000;C=2000;Period=10000;
                Count=indserangan;
                wave2=wave;
                for ii=1:length(wave)
                    z=mod(ii,Period+Count);
                    if z>=0 && z<Count
                        wave2(ii)=wave(z+DI);
                        
                    elseif z>=Count && z<DI
                        wave2(ii)=wave(z);
                    elseif z>=DI && z<DI+Count
                        wave2(ii)=wave(z-DI);
                    elseif z>=DI+C && z<Period
                        wave2(ii)=wave(z);
                    end
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==43 %CutSample
                %                 DI=6000;C=2000;Period=10000;
                RemoveNumber=7; Remove=1000;
                %                 Count=indserangan(2);
                
                wave2=zeros(size(wave));
                for ii=1:length(wave)
                    z=mod(ii,Remove-RemoveNumber);
                    wave2(ii)=wave(z+RemoveNumber);
                    
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
                
            elseif jenis==44 %ZeroCross
                ZeroCross=1000;
                
                wave2=wave;
                for ii=1:length(wave)
                    if abs(wave(ii))<ZeroCross/2^15
                        wave2(ii)=0;
                        
                    end
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==45 %ZeroLength1 Untuk kasus stereo, masing2 channel dinolkan tergantung masing2 channel
                Z=10;
                wave2=wave;
                for ii=1:length(wave)
                    if wave(ii)==0
                        wave2([ii:ii+Z])=0;
                        
                    end
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==46 %ZeroLength2 Untuk kasus stereo, semua channel dinolkan
                Z=10;
                wave2=wave;
                for ii=1:length(wave)
                    if wave(ii)==0
                        wave2([ii:ii+Z])=0;
                        
                    end
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==47 %ZeroRemove
                Z=10;
                %                 wave2=wave;
                for ii=1:length(wave)
                    if wave(ii)==0
                        if ii+1>length(wave)
                            break
                        end
                        wave2(ii)=wave(ii+1);
                    else
                        wave2(ii)=wave(ii);
                    end
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            elseif jenis==48 %Exchange
                
                wave2=wave;
                for ii=2:length(wave)-1
                    z=mod(ii,2);
                    if z==0
                        wave2(ii)=wave(ii+1);
                    elseif z==1
                        wave2(ii)=wave(ii-1);
                    end
                end
                    audiowrite(attackfile,wave2,Fs,'BitsPerSample',nbit)
                
            end
            
  end
end



function y=equalizer(x)
gaindb=[-6 6 -6 6 -6 6 -6 6 -6 6];
gain=10.^(gaindb/20);
fs=44100;
orde=1;
ft=2^5;
fc=ft/(fs/2);
[num1,den1]=butter(orde,fc);

ft=2^6;
fc=[ft-ft/2 ft+ft/2]/(fs/2);
[num2,den2]=butter(orde,fc);

ft=2^7;
fc=[ft-ft/2 ft+ft/2]/(fs/2);
[num3,den3]=butter(orde,fc);

ft=2^8;
fc=[ft-ft/2 ft+ft/2]/(fs/2);
[num4,den4]=butter(orde,fc);

ft=2^9;
fc=[ft-ft/2 ft+ft/2]/(fs/2);
[num5,den5]=butter(orde,fc);

ft=2^10;
fc=[ft-ft/2 ft+ft/2]/(fs/2);
[num6,den6]=butter(orde,fc);

ft=2^11;
fc=[ft-ft/2 ft+ft/2]/(fs/2);
[num7,den7]=butter(orde,fc);

ft=2^12;
fc=[ft-ft/2 ft+ft/2]/(fs/2);
[num8,den8]=butter(orde,fc);

ft=2^13;
fc=[ft-ft/2 ft+ft/2]/(fs/2);
[num9,den9]=butter(orde,fc);

ft=2^14;
fc=[ft]/(fs/2);
[num10,den10]=butter(orde,fc,'high');

y=zeros(size(x));
for i=1:10
    eval(['num' num2str(i) '=num' num2str(i) '*gain(i);']);
    eval(['y=y+filter(num' num2str(i) ',den' num2str(i) ',x);'])
    %     eval(['y=y+filtfilt(num' num2str(i) ',den' num2str(i) ',x);'])
    %     y1=filter(num1,den1,x);
end
y=y/max(abs(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Pink Noise Generation with MATLAB Implementation   %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov       07/30/13  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = pinknoise(N)

% function: y = pinknoise(N)
% N - number of samples to be returned in row vector
% y - row vector of pink (flicker) noise samples

% The function generates a sequence of pink (flicker) noise samples.
% Pink noise has equal energy in all octaves (or similar log bundles) of frequency.
% In terms of power at a constant bandwidth, pink noise falls off at 3 dB per octave.

% difine the length of the vector
% ensure that the M is even
if rem(N,2)
    M = N+1;
else
    M = N;
end

% generate white noise
x = randn(1, M);

% FFT
X = fft(x);

% prepare a vector for 1/f multiplication
NumUniquePts = M/2 + 1;
n = 1:NumUniquePts;
n = sqrt(n);

% multiplicate the left half of the spectrum so the power spectral density
% is proportional to the frequency by factor 1/f, i.e. the
% amplitudes are proportional to 1/sqrt(f)
X(1:NumUniquePts) = X(1:NumUniquePts)./n;

% prepare a right half of the spectrum - a copy of the left one,
% except the DC component and Nyquist frequency - they are unique
X(NumUniquePts+1:M) = real(X(M/2:-1:2)) -1i*imag(X(M/2:-1:2));

% IFFT
y = ifft(X);

% prepare output vector y
y = real(y(1, 1:N));

% ensure unity standard deviation and zero mean value
y = y - mean(y);
yrms = sqrt(mean(y.^2));
y = y/yrms;

% end



