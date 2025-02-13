% Function to take data from Synapse software (TDT) and change to matlab

function outputMat = saveFPmat(dateId, mouseId, synapseId);

tankId = [mouseId, '-', dateId, '-',synapseId];

cd(['folder location', dateId]); %CHANGE to folder with synapse file
temp = TDTbin2mat(['folder location\', dateId,'\', tankId]);

Ca = temp.streams.Dv1S.data; % calcium signal
input = temp.streams.Trig.data; % TTL signal
iso = temp.streams.Dv2S.data; % isobestic singal
Ch4 = temp.streams.Fi2r.data; % floating channel

% create plot of all traces
figure; 
subplot(3,1,1); plot(Ca);
subplot(3,1,2); plot(input);
subplot(3,1,3); plot(iso);
% figure; plot(Ch4)

filelocation = ['file location\', dateId, '\', mouseId,'-', dateId,'-', synapseId]; % CHANGE
cd(filelocation);

clearvars -except Ca iso input tankId ;
save([tankId]);

outputMat.Ca = Ca;
outputMat.input = input;
outputMat.iso = iso;
outputMat = outputMat;

end


