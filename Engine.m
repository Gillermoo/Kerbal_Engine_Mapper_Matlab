classdef Engine
    %Engine a rocket engine used by the KSP functions
    
    properties
        name string
        mass double
        Tasl double
        Tv double
        ISPasl double
        ISPv double
        cost double
        isRadial logical
        fuelType string
        builtInFuel double
        
    end
    
    methods (Static)
        function allowEngines = setupEngines(whichEngines)
            %mask = cell2mat(whichEngines(:,2))==0;
            %whichEngines(mask,:) = [];
            %setupEngines generates an array of engines based on which ones
            %are set to 1 in which engines
            [numE,~,rawE] = xlsread('RFEngines.csv');
            [numSRB,~,rawSRB] = xlsread('SRB.csv');
            allowEngines = [];
            for i = 1:length(numE(:,1))
                thisName = rawE{i,1};
                idx = find(strcmpi(whichEngines(:,1),thisName));
                if ~isempty(idx) && whichEngines{idx,2}==1
                    %generate a non-srb fueled rocket engine
                    isRadial = strcmpi('Radial mounted',rawE{i,2});
                    
                    if strcmpi('Nerv',rawE{i,1})
                        type = 'LF';
                    elseif strcmpi('Dawn',rawE{i,1})
                        type = 'Xenon';
                    else
                        type = 'LFOX';
                    end
                    
                    if strcmpi('Twin-Boar',rawE{i,1})
                        BIF = 6400;
                    else
                        BIF = 0;
                    end
                    allowEngines = [allowEngines,Engine(rawE{idx,1},numE(idx,2),numE(idx,6),...
                       numE(idx,7), numE(idx,11), numE(idx,12), numE(idx,1),...
                       isRadial,type,BIF)]; 
                end
            end
            
            for i = 1:length(numSRB(:,1))
                thisName = rawSRB{i,1};
                idx = find(strcmpi(whichEngines(:,1),thisName))-length(numE(:,1));
                if ~isempty(idx) && whichEngines{idx+length(numE(:,1)),2}==1
                    %generate a non-srb fueled rocket engine
                    isRadial = strcmpi('Radial mounted',rawSRB{i,2});
                    BIF = numSRB(idx,7);
                    allowEngines = [allowEngines,Engine(rawSRB{idx,1},numSRB(idx,3),numSRB(idx,8),...
                       numSRB(idx,9), numSRB(idx,14), numSRB(idx,15),numSRB(idx,1),...
                       isRadial,'SRB',BIF)]; 
                end
            end
        end
        
    end
    
    methods
        function obj = Engine(name,mass,Tasl,Tv,ISPasl,ISPv,cost,isRadial,fuelType,builtInFuel)
            obj.name = name;
            obj.mass = mass;
            obj.Tasl = Tasl;
            obj.Tv = Tv;
            obj.ISPasl = ISPasl;
            obj.ISPv = ISPv;
            obj.cost = cost;
            obj.isRadial = isRadial;
            obj.fuelType = fuelType;
            obj.builtInFuel = builtInFuel;
        end
    end
end

