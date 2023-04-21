function [paretoData] = paretoPoints(data,goals,varargin)
%pareto determines the pareto optimal points for the rows in a dataset
if isempty(varargin)
    indicies = 1:length(data(:,1));
else
    indicies = varargin{1};
end
if length(data(:,1))>1
    [~,bestidx] = max(data.*goals,[],1);
    bestidx = unique(bestidx);
    paretoData = [data(bestidx,:),indicies(bestidx)'];
    if length(bestidx)~=1
        paretoData = [data(bestidx,:),indicies(bestidx)'];
        data(bestidx,:) = [];
        indicies(bestidx) = [];
        mask = zeros([length(data(:,1)),1]);
        adjustedData = data.*goals;
        for i = 1:length(goals)
            idx = 1:length(goals);
            idx(i) = [];
            tempMask = any(paretoData(i,idx).*goals(:,idx)>adjustedData(:,idx),2);
            mask = mask | tempMask;
        end
        data(mask,:) = [];
        indicies(mask) = [];
        if ~isempty(data)
            paretoData = [paretoData;paretoPoints(data,goals,indicies)];
        end
    end
    
    
else
    paretoData = [data,indicies];
end
end

