clear all
%% OFDM with MIMO Simulation
NUM=10;
PAPR=zeros(1,NUM);
ii=[1,4,9]; %setting k 
nn=size(ii,2);
RE=cell(1,nn+1);
M=cell(1,nn+1);
iteration=20;
iteration2=20;
cc=0;
cc_N=[];
cc_A=[];
for i=1:nn+1
    M{i}=zeros(1,iteration+1);
end
count=1;
flag=1;
% Create a QPSK modulator and demodulator pair.
qpskMod = comm.QPSKModulator;
qpskDemod = comm.QPSKDemodulator;
NU=0;
%%
ofdmMod = comm.OFDMModulator('FFTLength',28,'PilotInputPort',true,...
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
%  60 RBs
nframes = 60;
while(1)
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
    dataOFDM=reshape(dataOFDM,12,2);
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
w=ones(size(A,2),1);
alpha=sum((A*w).*conj(A*w))/(size(A*w,1));
PAPR(NU+1)=max((A*w).*conj(A*w))/alpha;
if PAPR(NU+1)>12.75
    continue
end

w=cell(1,nn+3);


for i=1:nn
[w{i},S2,m]=NearestVector_new(A,iteration,alpha,ii(i));
average=sum(S2.*conj(S2))/size(S2,1);
RE{i}(NU+1)=max(S2.*conj(S2))/average;
M{i}=M{i}+m;
end

[w{nn+1},m]=orimethod(A,iteration2);
M{nn+1}=M{nn+1}+m;
S2=A*w{nn+1};
average=sum(S2.*conj(S2))/size(S2,1);
RE{nn+1}(NU+1)=max(S2.*conj(S2))/average;
NU=NU+1;
fprintf("iteration: %d \n",NU);
if NU==NUM
    break
end

end


for i=1:nn+1
    M{i}=M{i}./NUM;
end

for i=1:nn
    RE{i}=sort(RE{i});
end

mark=[1:9000/3:9000,9001:900/3:9900,9901:90/3:9990,9991:9/3:10000];
mark4=[4000:9000/3:9000,9501:900/3:9900,9951:90/3:9990,9995:9/3:10000];
mark2=1:4:21;
mark3=2:4:21;

figure(1)
G=semilogy(sort(PAPR),1-linspace(0,1,NUM),'k-');hold on

H(1)=semilogy(sort(RE{1}),1-linspace(0,1,NUM),'ro--','MarkerIndices',mark);hold on
H(2)=semilogy(sort(RE{2}),1-linspace(0,1,NUM),'mx--','MarkerIndices',mark);hold on
H(3)=semilogy(sort(RE{3}),1-linspace(0,1,NUM),'bd--','MarkerIndices',mark4);hold on
H(4)=semilogy(sort(RE{4}),1-linspace(0,1,NUM),'--');hold on


grid on
xlabel('PAPR(dB)');
ylabel('CCDF');
legend([G,H(1:4)],'original','k=1','k=2','k=3','UC-CMA')

figure(2)
H1=plot(0:iteration,M{1},'ro--','MarkerIndices',mark2);hold on;
H2=plot(0:iteration,M{2},'gs--','MarkerIndices',mark2);hold on;
H3=plot(0:iteration,M{3},'bd--','MarkerIndices',mark3);hold on;
H4=plot(0:iteration,M{4},'--');hold on;

legend([H1,H2,H3,H4],'k=1','k=2','k=3','UC-CMA');
xlabel("iterations")
ylabel("PAPR")

%%

function [w,M]=orimethod(A,iteration)
tstart=tic;
N=size(A,2);
w=ones(N,1);
alpha=sum((A*w).*conj(A*w))/(size(A*w,1));
tp=0.005/(alpha)^2;
ss=(A*w).*conj(A*w);
M=10*log10(max(ss)/mean(ss));
w=ones(N,1);

for i=1:iteration
    s=A*w;
    e=(s.*conj(s))-alpha;
    se=s.*e;
    w=w-tp*A'*se;
    w=w./abs(w);
    M=[M,10*log10(max((A*w).*conj(A*w))/alpha)];

end
tend=toc(tstart);
fprintf("Elapsed time is %f seconds with original method \n",tend)
hold on
end




function AB = kr(A,B)
%% KR Khatri-Rao product

[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=B(:,f)*A(:,f).';
   AB(:,f)=ab(:);
end
end
