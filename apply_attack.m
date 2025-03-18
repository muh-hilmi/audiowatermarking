function xwn = apply_attack(xw, fs, serangan, nbit)
    allattack_audio_stirmark(serangan, 'temphost.wav', 'temp_attack.wav', nbit); % Placeholder
    [xwn, ~] = audioread('temp_attack.wav');
end