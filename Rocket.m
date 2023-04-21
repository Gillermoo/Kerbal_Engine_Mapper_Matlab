classdef Rocket
    %Rocket KSP Rocket
    
    properties
        stages RocketStage
        cost double
        mass double
        dv double
    end
    
    methods
        
        function obj = Rocket(varargin)
            %Rocket generates a rocekt either based on an array of stages,
            %or a set of inputs to find an optimal design
            if isa(varargin{1}, 'RocketStage')
                %generate rocket based on specified stages
                obj.stages = varargin{1};
                obj.cost = 0;
                obj.mass = 0;
                obj.dv = 0;
                for i = 1:size(obj.stages,1)
                    obj.cost = obj.cost + obj.stages(i).cost;
                    obj.mass = obj.mass + obj.stages(i).mass;
                    obj.dv = obj.dv + obj.stages(i).dv;
                end
            else
                %generate optimal rocket stages base based on input
                %requirements
                stages = [];
                inData = varargin{1};
                %varargin{2} is allowed engines
                %varargin{3} is fuel tanks data
                cost = 0;
                mass = 0;
                dv = 0;
                
                for i = 1:length(inData)
                    data = inData{i};
                    stageInfo = struct;
                    while ~isempty(data)
                        [str1,data] = strtok(data,' :');
                        [str2,data] = strtok(data,' :');
                        if strcmpi(str1,'t') || strcmpi(str1,'cond')
                            stageInfo.(str1) = str2;
                        else
                            stageInfo.(str1) = str2num(str2);
                        end
                        
                    end
                    if isempty(stages)
                        %if not stage has been run yet, add the first stage
                        %to the top of the stack
                        switch lower(stageInfo.t)
                            case 'pl'
                                stages = RocketStage.PLStage(stageInfo.pl);

                        end
                    else
                        newstages = [];
                        totalEngines = length(varargin{2}).*(length(varargin{2})-1);
                        for j = 1:length(stages(1,:))
                            switch lower(stageInfo.t)
                                case 'pl'
                                    newstages = [newstages,RocketStage.PLStage(stageInfo.pl)];
                                case 'lin'
                                    newstages = [newstages,RocketStage.linearStage(stageInfo,stages(:,j),varargin{2},varargin{3})];

                                case 'asp'
                                    newstages = [newstages,RocketStage.aspStage(stageInfo,stages(:,j),varargin{2},varargin{3})];
                            end
                        end
                        
                        if ~isempty(newstages)
                            massArray = zeros([1,length(newstages(1,:))]);
                            costArray = zeros([1,length(newstages(1,:))]);
                            for j = 1:length(newstages(1,:))
                                thismass = 0;
                                thiscost = 0;
                                for k = 1:length(newstages(:,1))
                                    thismass = thismass + newstages(k,j).mass;
                                    thiscost = thiscost + newstages(k,j).cost;
                                end
                                %disp(stageInfo)
                                massArray(j) = thismass;
                                costArray(j) = thiscost;
                            end
                            %scatter(massArray,costArray)

                            %if nextstage has no engines, keep all rockets
                            if i == length(inData)
                                
                                paretoData = paretoPoints([massArray',costArray'],[-1,-1]);
                                %scatter(paretoData(:,1),paretoData(:,2));
                                stages = newstages(:,paretoData(:,end));
                            else
                                test = inData{i+1};
                                if contains(test, 'center:0')
                                    stages = newstages;
                                else
                                    paretoData = paretoPoints([massArray',costArray'],[-1,-1]);
                                    stages = newstages(:,paretoData(:,end));
                                end
                            end
                                
                        else
                            error('It appears no engines can meet your requirements');
                        end
                    end
                    
                    %cost = cost + stages(i).cost;
                    %mass = mass + stages(i).mass;
                    %dv = dv + stages(i).dv;
                end
                %obj.stages = stages;
                obj = Rocket.empty;
                for i = 1:size(stages,2)
                    newobj = Rocket(stages(:,i));
                    obj = [obj,newobj];
                end
            end
            %scatter(costArray,massArray,15,'filled')
            
        end
        
        function [tbl, data, names, rows] = toString(obj)
            names = {'Mass','Delta_V','Cost','Engine','Type','Liquid_Fuel','Xenon'};
            idx = 1;
            data = {};
            rows = {};
            numStages = length(obj.stages);
            for i = 1:length(obj.stages)
                s = obj.stages(i);
                data{idx,1} = s.mass;
                if(isempty(s.engines))
                    data{idx,4} = [];
                else
                    data{idx,4} = s.engines.name;
                end
                
                data{idx,2} = s.dv;
                data{idx,5} = s.type;
                data{idx,6} = s.lFuel;
                data{idx,7} = s.xenon;
                data{idx,3} = s.cost;
                if strcmp(s.type,'Payload')
                    rows{idx} = sprintf('%0.0f Payload', numStages-idx);
                else
                    rows{idx} = sprintf('%0.0f Stage', numStages-idx);
                end
                
                idx = idx+1;
            end
            rows{idx} = sprintf('Total', numStages-idx);
            data{idx,1} = obj.mass;
            data{idx,2} = obj.dv;
            data{idx,3} = obj.cost;
            
            outstr = cell2table(data,'VariableNames',names,'RowNames',rows);
            %varargout = {data,'VariableNames',names,'RowNames',rows};
            tbl = cell2table(data,'VariableNames',names,'RowNames',rows);
            for i = 1:length(data(:,1))
                for j = 1:length(data(1,:))
                    if(isnumeric(data{i,j}))
                        data{i,j} = num2str(data{i,j});
                    else
                        data{i,j} = char(data{i,j});
                    end
                end
            end
            
        end
        
    end
   
end

