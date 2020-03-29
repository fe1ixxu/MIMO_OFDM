clear all
%% OFDM with MIMO Simulation
NUM=1;
PAPR=zeros(1,NUM);
RE=zeros(1,NUM);
for NU=1:NUM

% Create a QPSK modulator and demodulator pair.
qpskMod = comm.QPSKModulator;
qpskDemod = comm.QPSKDemodulator;
%%
ofdmMod = comm.OFDMModulator('FFTLength',32,'PilotInputPort',true,...
    'PilotCarrierIndices',cat(3,[8;22],...
    [9;23]),'InsertDCNull',true,...
    'NumTransmitAntennas',2);

ofdmDemod = comm.OFDMDemodulator(ofdmMod);
ofdmDemod.NumReceiveAntennas = 2;

% showResourceMapping(ofdmMod)
ofdmModDim = info(ofdmMod);

numData = ofdmModDim.DataInputSize(1);   % Number of data subcarriers
numSym = ofdmModDim.DataInputSize(2);    % Number of OFDM symbols
numTxAnt = ofdmModDim.DataInputSize(3);  % Number of transmit antennas
%%
%  10 RBs
nframes = 32*4;
data = randi([0 3],nframes*numData,numSym,numTxAnt);
%%
% Apply QPSK modulation
modData = qpskMod(data(:));
modData = reshape(modData,nframes*numData,numSym,numTxAnt);
%%
% Create an error rate counter.
errorRate = comm.ErrorRate;
%%
% Simulate the OFDM system 
A2=[];
B=[];
W=[];
U_svd=struct([]);
for k = 1:nframes

    % Find row indices for kth OFDM frame
    indData = (k-1)*ofdmModDim.DataInputSize(1)+1:k*numData;

    % Generate random OFDM pilot symbols
    pilotData = complex(rand(ofdmModDim.PilotInputSize), ...
        rand(ofdmModDim.PilotInputSize));
    pilotData =reshape(pilotData,2,2);
    % Modulate QPSK symbols using OFDM
    dataOFDM =modData(indData,:,:);
    dataOFDM=reshape(dataOFDM,16,2);
%     dataOFDM=[dataOFDM(1:7,:);pilotData;dataOFDM(8:end,:)];
    B=blkdiag(B,dataOFDM.');
    B1=B;
    chGain=eye(2); 
    
    [ls,sv,rs]=svd(chGain);
     W=[W,rs]; %precoding
     U_svd(k).V=rs;
     U_svd(k).ch=chGain;
     U_svd(k).U=ls; %sovle precoding
end

TotalSize=numData;
N=TotalSize*nframes;
F=ifft(eye(N));
B=B*F';
A=kr(B.',W);
iteration=10;
w=ones(size(A,2),1);
alpha=sum((A*w).*conj(A*w))/(size(A*w,1));
PAPR(NU)=max((A*w).*conj(A*w))/alpha;

% testing num
ii=[1,4,9];
nn=length(ii);
w=cell(1,nn+1);
ww=cell(1,nn+1);
S=cell(1,nn+1);

for i=1:nn
[w{i},S2,m]=NearestVector_new(A,iteration,alpha,ii(i));
ww{i}=diag(w{i});
S{i}=ww{i}*B*fft(eye(N))';
end
w{nn+1}=orimethod(A,iteration);
ww{nn+1}=diag(w{nn+1});
S{nn+1}=ww{nn+1}*B*fft(eye(N))';

%%
errorrate=cell(1,nn+1);
for i=1:nn+1
vecs=A*w{i};
st=1;
for numBER=1:st:15
    
for k = 1:nframes
    indData = (k-1)*ofdmModDim.DataInputSize(1)+1:k*numData;
    indData2 = 2*(k-1)*TotalSize+1:2*k*TotalSize;
    
    % Apply least squares solution to remove effects of fading channel
    

    % Demodulate OFDM data
    
    SigTime=S{i}(2*k-1:2*k,indData);
    wwtemp=ww{i}(2*k-1:2*k,2*k-1:2*k);
%     SigTime=fft(SigTime.').';
    % channel
    SigTime=U_svd(k).U'*U_svd(k).ch*SigTime;
    SigTime=awgn(SigTime,numBER);
    SigTime=wwtemp\SigTime;
    SigTime=SigTime.';
    % Demodulate QPSK data
    receivedData = qpskDemod(SigTime(:));
    
    % Compute error statistics
    dataTmp = data(indData,:,:);
    errors = errorRate(dataTmp(:),receivedData);
end
errorrate{i}=[errorrate{i} errors(1)];
reset(errorRate);
end
end


end



H(1)=semilogy(1:st:numBER,errorrate{1},'ro-');hold on;
H(2)=semilogy(1:st:numBER,errorrate{2},'mx-');hold on;
H(3)=semilogy(1:st:numBER,errorrate{3},'bd-');hold on;
H(4)=semilogy(1:st:numBER,errorrate{4},'-');hold on;

ylim([10^-3 1]);
xlabel('SNR(dB)')
ylabel('BER')
grid on
legend(H(1:4),'k=1','k=2','k=3','UC-CMA')

%%

function w=orimethod(A,iteration)

N=size(A,2);
w=ones(N,1);
alpha=sum((A*w).*conj(A*w))/(size(A*w,1));
tp=0.05/(alpha)^2;
M=[];
w=ones(N,1);

for i=1:iteration
    M=[M,10*log10(max((A*w).*conj(A*w))/alpha)];
    s=A*w;
    e=(s.*conj(s))-alpha;
    se=s.*e;
    w=w-tp*A'*se;
    w=w./abs(w);

end

% plot(1:iteration,M,'b')
hold on
end




function AB = kr(A,B)
%% KR Khatri-Rao product

[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' Error in kr.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=B(:,f)*A(:,f).';
   AB(:,f)=ab(:);
end
end
