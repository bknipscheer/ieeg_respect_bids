function samplocs = findtrigger(data, stimnum,fs, sampstart,sampend)

if median(data(stimnum(1),sampstart:sampend)) > median(data(stimnum(2),sampstart:sampend))
    anode = stimnum(2);% with positive saturation right after stimulation
    cathode = stimnum(1);% with negative saturation right after stimulation
else
    anode = stimnum(1);% with positive saturation right after stimulation
    cathode = stimnum(2);% with negative saturation right after stimulation
end

dt_data(1,:) = diff(data(anode,sampstart:sampend));
dt_data(2,:) = diff(data(cathode,sampstart:sampend));
SD = std(dt_data,[],2);

locsdtpos = []; locsdtneg = [];
if any(dt_data(1,:)>10*SD(1))
    [~,locsdtpos] = findpeaks(dt_data(1,:),'MinPeakDistance',fs,'MinPeakHeight',10*SD(1));
end
if any(-1*dt_data(1,:)>10*SD(1))
    [~,locsdtneg] = findpeaks(-1*dt_data(1,:),'MinPeakDistance',fs,'MinPeakHeight',10*SD(1));
end

locsdt2pos = []; locsdt2neg = [];
if any(dt_data(2,:)>10*SD(2))
    [~,locsdt2pos] = findpeaks(dt_data(2,:),'MinPeakDistance',fs,'MinPeakHeight',10*SD(2));
end
if any(-1*dt_data(2,:)>10*SD(2))
    [~,locsdt2neg] = findpeaks(-1*dt_data(2,:),'MinPeakDistance',fs,'MinPeakHeight',10*SD(2));
end

locsall = sort([locsdtpos locsdtneg, locsdt2pos, locsdt2neg]);
distances = pdist2(locsall',locsall');

within1 = distances < round(0.5*fs); % group all locs that are within 500ms
numelements = numel(locsall);
% use only lower triangle part
t = logical(triu(ones(numelements,numelements),1)); 
within1(t) = 0;
[labeledGroups,numGroups ] = bwlabel(within1,4);

for k = 1 : numGroups
  [rows, columns] = find(labeledGroups == k);
  indexes = unique([rows, columns]);
  if k>1
      if locsall(indexes(1))-samplocs(k-1) < round(6*fs) % WARNING: assumes that stimuli are applied maximally every 5s
          samplocs(k) = min(locsall(indexes));
      else
          break
      end
  else
      samplocs(k) = min(locsall(indexes));
  end
end
