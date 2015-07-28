% N = 8;  % Length of Walsh (Hadamard) functions

HadIdx = 0:N-1;                          % Hadamard index
M = log2(N)+1;                           % Number of bits to represent the index
 
binHadIdx = fliplr(dec2bin(HadIdx,M))-'0'; % Bit reversing of the binary index
binSeqIdx = zeros(N,M-1);                  % Pre-allocate memory
for k = M:-1:2
    % Binary sequency index
    binSeqIdx(:,k) = xor(binHadIdx(:,k),binHadIdx(:,k-1));
end
SeqIdx = binSeqIdx*pow2((M-1:-1:0)');    % Binary to integer sequency index

scene_hp = scene_hp(:,SeqIdx+1);
scene_hn = scene_hn(:,SeqIdx+1);
sensor_hp = sensor_hp(:,SeqIdx+1);
sensor_hn = sensor_hn(:,SeqIdx+1);

%% Test
% hadamardMatrix = hadamard(N);
% fwht_hadamard = fwht(eye(N)*N,N,'hadamard');
% [hadamardMatrix-fwht_hadamard]
% 
% walshMatrix = hadamardMatrix(SeqIdx+1,:); % 1-based indexing
% fwht_walsh = fwht(eye(N)*N,N);
% [walshMatrix-fwht_walsh]