classdef RocketStage
    %RocketStage KSP Rocket Stage
    
    properties
        mass double
        engines Engine
        numEngines uint8
        cost double
        type string
        dv double
        fuel double
        lFuel double
        ox double
        xenon double
    end
    
    methods
        
        function obj = RocketStage()
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.mass = 100;
            %obj.engines = '2 x small one';
            obj.numEngines = 1;
            obj.cost = 69;
            obj.type = 'linear';
            obj.dv = 2000.12431231312313;
            obj.fuel = 123;
            obj.lFuel = 100;
            obj.ox = 23;
        end
        
        function [] = toString(obj,varargin)
            %toString Prints out the stage to a string
            if varargin{1}==1
                str1 = sprintf('| %15s | %15s | %15s | %15s | %15s |','Mass','Engine','Delta-V','Type','Liquid Fuel');
                disp(str1)
            end
            str2 = sprintf('| %15.0f | %15s | %15.0f | %15s | %15.0f |',obj.mass,obj.engines,obj.dv,obj.type,obj.lFuel);
            disp(str2);
        end
        
    end
    
    methods (Static)
        
        function obj = PLStage(pl)
            %PLStage generates a payload stage with not fuel or engines
            obj = RocketStage;
            obj.mass = pl;
            obj.numEngines = 0;
            obj.cost = 0;
            obj.type = 'Payload';
            obj.dv = 0;
            obj.fuel = 0;
            obj.lFuel = 0;
            obj.ox = 0;
            obj.xenon = 0;
        end
        
        function [h] = plotLinearStage(stageInfo,mass0,engines,tanks)
            g = 9.8;
            Ps = .11252; %percent of fueltank that is structural
            switch stageInfo.cond
                
                case 'asl'
                    cond = 1;
                case 'v'
                    cond  = 2;
            end
            outStages = [];
            for eng = engines
                %loop through the number of engines
                %disp(eng.name)
                T = [eng.Tasl,eng.Tv];
                ISP = [eng.ISPasl,eng.ISPv];
                %                  LF, LFOX, Xenon, SRB
                FuelDensities = [.005, .005, .001, .0075];
                
                switch eng.fuelType
                    case 'LF'
                        densityIdx = 1;
                        totfuel = tanks.LF(:,end);
                        tankCost = tanks.LF(:,1);
                    case 'LFOX'
                        densityIdx = 2;
                        totfuel = tanks.RF(:,end) + tanks.RF(:,end-1);
                        tankCost = tanks.RF(:,1);
                    case 'Xenon'
                        densityIdx = 3;
                        totfuel = tanks.Xenon(:,end);
                        tankCost = tanks.Xenon(:,1);
                    case 'SRB'
                        densityIdx = 4;
                end
                
                for numEng = 1:stageInfo.e
                    %skip iteration if a single radial engine is chosen
                    if eng.isRadial && numEng==1
                        continue
                    end
                    
                    
                    if densityIdx == 4
                        Mff = exp(stageInfo.dv./(ISP(cond).*g));
                        M0 = Mff.*(numEng.*eng.mass+mass0);
                        Mf = M0-mass0-numEng.*eng.mass;
                        fuel = Mf/FuelDensities(densityIdx);
                        fuel(fuel>numEng.*eng.builtInFuel) = nan;
                        Mt = numEng.*eng.mass+fuel.*FuelDensities(densityIdx);
                        obj = RocketStage;
                        obj.mass = numEng.*eng.mass + Mt;
                        obj.engines = eng;
                        obj.numEngines = numEng;
                    else
                        %add fuel
                        %mass Fuel Fraction
                        Mff = exp(stageInfo.dv./(ISP(cond).*g));
                        %Mass of fuel tanks at start
                        Mt = -(numEng.*eng.mass+mass0).*(Mff-1)./(Mff.*Ps-1);
                        %if engine cant meet TWR skip
                        Mf = (1-Ps)*Mt;
                        fuel = Mf/FuelDensities(densityIdx);
                        if strcmp(eng.name,'Twin-Boar')
                            if fuel>=6400*numEng
                                fuel = fuel-6400*numEng;
                            else
                                continue
                            end

                        end
                        obj = RocketStage;
                        obj.mass = Mt+numEng*eng.mass;
                        obj.engines = eng;
                        obj.numEngines = numEng;
                    end
                    %check to see if the solid fuel booster has enough fuel
                    %for the required Delta-V
                    %Make sure twin boar is using up all its fuel
                    
                    
                    
                    
                    if densityIdx == 4
                        obj.cost = numEng*eng.cost.*ones(size(Mf));
                        fuel(fuel>eng.builtInFuel) = nan;
                    else
                        %disp(eng)
                        
                            
                        cpp = tankCost./totfuel;
                        cppIDX = knnsearch(totfuel,fuel(:));
                        CostPerFuel = cpp(cppIDX);
                        %CostPerFuel = interp1(totfuel,tankCost./totfuel,1000, 'nearest','extrap');
                        CostPerFuel = reshape(CostPerFuel, size(fuel));
                        obj.cost = numEng.*eng.cost + CostPerFuel.*fuel;
                        obj.cost(fuel<0) = nan;
                    end
                    obj.cost(fuel<0) = nan;
                    obj.cost(isnan(fuel)) = nan;
                    obj.mass(fuel<0) = nan;
                    obj.mass(isnan(fuel)) = nan;
                    obj.cost(obj.mass<0) = nan;
                    obj.cost(numEng.*T(cond)<stageInfo.twr.*(mass0 + numEng.*eng.mass+Mt)*g) = nan;
                    obj.mass(numEng.*T(cond)<stageInfo.twr.*(mass0 + numEng.*eng.mass+Mt)*g) = nan;
                    obj.type = 'Linear';
                    obj.dv = stageInfo.dv;
                    obj.fuel = fuel;
                    
                    switch eng.fuelType
                        
                        case 'LF'
                            obj.lFuel = fuel;
                            obj.ox = 0;
                            obj.xenon = 0;
                        case 'LFOX'
                            obj.lFuel = fuel*.45;
                            obj.ox = fuel*.55;
                            obj.xenon = 0;
                        case 'Xenon'
                            obj.xenon = fuel;
                            obj.lFuel = 0;
                            obj.ox = 0;
                        case 'SRB'
                            obj.xenon = 0;
                            obj.lFuel = 0;
                            obj.ox = 0;
                    end
                    obj.mass(obj.mass<0) = nan;
                    if densityIdx == 4
                        unitsFuel = Mf.* FuelDensities(4);
                        obj.mass(unitsFuel>obj.engines.builtInFuel) = nan;
                        obj.cost(unitsFuel>obj.engines.builtInFuel) = nan;
                    end
                    outStages = [outStages,obj];
                    %if engine meets TWR then dont add more engines
                end
            end
            h = figure('units','normalized','outerposition',[0 0 1 1]);
            [ha, pos] = tight_subplot(1, 2, [.01 .03],[.1 .05],[.05 .01]);
            for type = 1:2
                if type == 1
                    bestval = 1e20.*ones(size(Mf));
                    bestmass = 1e20.*ones(size(Mf));
                    bestidx = ones(size(Mf)).*nan;
                    for i = 1:length(outStages)
                        bestidx(outStages(i).cost < bestval) = i;
                        bestmass(outStages(i).cost < bestval) = outStages(i).mass(outStages(i).cost < bestval);
                        bestval(outStages(i).cost < bestval) = outStages(i).cost(outStages(i).cost < bestval);
                        
                    end
                else
                    bestval = 1e20.*ones(size(Mf));
                    bestidx = ones(size(Mf)).*nan;
                    for i = 1:length(outStages)
                        bestidx(outStages(i).mass < bestval) = i;
                        bestmass(outStages(i).mass < bestval) = outStages(i).mass(outStages(i).mass < bestval);
                        bestval(outStages(i).mass < bestval) = outStages(i).mass(outStages(i).mass < bestval);
                        
                    end
                end
                
                
                axes(ha(type));
                %semilogy(1,1); hold on
                map = colorcube(length(outStages)+2);
                sumMap = sum(map,2);
                map(sumMap==0,:) = [];
                sumMap = sum(map,2);
                map(sumMap==3,:) = [];
                %[C,h] = contourf(obj.dv,mass0,bestidx);
                %set(h,'LineColor','none')
                bestidx(bestidx==0) = nan;
                imagesc(flipud(bestidx));hold on;
                colormap(map);
                caxis([1,length(outStages)]);
                grid on
                
                spacing = 10^floor(log10(max(max(obj.dv))-min(min(obj.dv))));
                dvs = 0:spacing:max(max(obj.dv));
                dvidx = interp1(obj.dv(1,:),1:length(obj.dv(1,:)),dvs);
                dvs(isnan(dvidx)) = [];
                dvidx(isnan(dvidx)) = [];
                xticks(dvidx);
                xticklabels(dvs);

                exponent = log10(mass0(:,1));
                mine = floor(min(exponent));
                maxe = ceil(max(exponent));
                
                m = logspace(mine,maxe,maxe-mine+1);

                midx = interp1(mass0(:,1),fliplr(1:length(mass0(:,1))),m);
                mm = m.*[1;2;3;4;5;6;7;8;9];
                m(isnan(midx)) = [];
                midx(isnan(midx)) = [];
                yticks(fliplr(midx));
                yticklabels(fliplr(m));

                set(gca,'YMinorTick','on')
                ax = gca;
                
                mm = mm(:);
                mmidx = interp1(mass0(:,1),1:length(mass0(:,1)),mm);
                mm(isnan(mmidx)) = [];
                mmidx(isnan(mmidx)) = [];
                ax.YAxis.MinorTickValues = flipud(length(mass0(:,1))-mmidx);
                if(length(m)<2)
                    mine = min(exponent);
                    maxe = max(exponent);
                    m = logspace(mine,maxe,ceil(maxe-mine)+1);
                    midx = interp1(mass0(:,1),fliplr(1:length(mass0(:,1))),m);
                    m(isnan(midx)) = [];
                    midx(isnan(midx)) = [];
                    yticks(fliplr(midx));
                    yticklabels(fliplr(m));
                end
                
                spacing = 10^floor(log10(max(max(obj.dv))-min(min(obj.dv))))/5;
                dvs = 0:spacing:max(max(obj.dv));
                dvidx = interp1(obj.dv(1,:),1:length(obj.dv(1,:)),dvs);
                dvs(isnan(dvidx)) = [];
                dvidx(isnan(dvidx)) = [];
                ax.XAxis.MinorTickValues = dvidx;
                grid minor
                grid on
                xlabel('Delta-V (m/s)');
                ylabel('Payload Mass (tons)');
                marked = zeros(size(Mf));
                for i = 1:length(Mf(:,1))
                    for j = 1:length(Mf(1,:))
                        if(marked(i,j) == 0)
                            if(isnan(bestidx(i,j)))
                                continue;
                            end
                            allOfThis = find(bestidx == bestidx(i,j));
                            
                            dvs = obj.dv(allOfThis);
                            masses = mass0(allOfThis);
                            
                            centerDv = sum(dvs)./length(dvs);
                            centerMass = sum(log10(masses))./length(masses);
                            
                            p = [centerDv,centerMass];
                            
                            cppIDX = knnsearch([dvs,log10(masses)],p);
                            
                            idx = cppIDX;
                            loc = allOfThis(idx);
                            m = mass0(loc);
                            dv = obj.dv(loc);
                            dvidx = interp1(obj.dv(1,:),1:length(obj.dv(1,:)),dv);
                            midx = length(mass0(:,1)) - interp1(mass0(:,1),1:length(mass0(:,1)),m);
                            if(sum(map(bestidx(i,j),:))< 1.25)
                                text(dvidx,midx,['\leftarrow ',char(num2str(outStages(bestidx(i,j)).numEngines)),' ', char(outStages(bestidx(i,j)).engines.name)],'FontWeight','bold','fontsize',10,'rotation', 30,'color','w');
                                scatter(dvidx,midx,30,'w');
                            else
                                text(dvidx,midx,['\leftarrow ',char(num2str(outStages(bestidx(i,j)).numEngines)),' ', char(outStages(bestidx(i,j)).engines.name)],'FontWeight','bold','fontsize',10,'rotation', 30);
                                scatter(dvidx,midx,30,'k');
                            end
                            scatter(dvidx,midx,30,map(bestidx(i,j),:),'filled');
                            marked(bestidx == bestidx(i,j)) = 1;
                        end
                    end
                end
                
                if(type == 1)
                    title(['Lowest Cost TWR: ',char(num2str(stageInfo.twr)),' condition: ',stageInfo.cond]);
                else
                    title(['Lowest Mass TWR: ',char(num2str(stageInfo.twr)),' condition: ',stageInfo.cond]);
                end
                
            end
            
            %saveas(gcf,['TWR ',char(num2str(stageInfo.twr)),' condition ',stageInfo.cond,' max engines ',char(num2str(stageInfo.e)),'.png']);
            
        end
        
        function outStages = linearStage(stageInfo,stages,engines,tanks)
            %LinearStage generates a rocket stage that is just a fuel tank
            %and engines with no stage separations
            g = 9.8;
            Ps = .11252; %percent of fueltank that is structural
            %get Mass0
            outStages = [];
            mass0 = 0;
            for i = 1:length(stages)
                mass0 = mass0 + stages(i).mass;
            end
            
            
            %loop through each engine
            switch stageInfo.cond
                
                case 'asl'
                    cond = 1;
                case 'v'
                    cond  = 2;
            end
            for eng = engines
                %loop through the number of engines
                %disp(eng.name)
                T = [eng.Tasl,eng.Tv];
                ISP = [eng.ISPasl,eng.ISPv];
                %                  LF, LFOX, Xenon, SRB
                FuelDensities = [.005, .005, .001, .0075];
                
                switch eng.fuelType
                    case 'LF'
                        densityIdx = 1;
                        totfuel = tanks.LF(:,end);
                        tankCost = tanks.LF(:,1);
                    case 'LFOX'
                        densityIdx = 2;
                        totfuel = tanks.RF(:,end) + tanks.RF(:,end-1);
                        tankCost = tanks.RF(:,1);
                    case 'Xenon'
                        densityIdx = 3;
                        totfuel = tanks.Xenon(:,end);
                        tankCost = tanks.Xenon(:,1);
                    case 'SRB'
                        densityIdx = 4;
                end
                
                for numEng = 1:stageInfo.e
                    %skip iteration if a single radial engine is chosen
                    if eng.isRadial && numEng==1
                        continue
                    end
                    %if engine cant meet TWR skip
                    if numEng*T(cond)<stageInfo.twr*(mass0 + numEng*eng.mass)*g
                        continue
                    end
                    %add fuel
                    %mass Fuel Fraction
                    Mff = exp(stageInfo.dv./(ISP(cond).*g));
                    %Mass of fuel tanks at start
                    Mt = -(numEng.*eng.mass+mass0).*(Mff-1)./(Mff.*Ps-1);
                    %if engine cant meet TWR skip
                    if Mt<0
                        continue
                    end
                    if numEng*T(cond)<stageInfo.twr*(mass0 + numEng*eng.mass+Mt)*g
                        continue
                    end
                    Mf = (1-Ps)*Mt;
                    fuel = Mf/FuelDensities(densityIdx);
                    %check to see if the solid fuel booster has enough fuel
                    %for the required Delta-V
                    if densityIdx == 4 && fuel> eng.builtInFuel
                        continue
                    end
                    %Make sure twin boar is using up all its fuel
                    if strcmp(eng.name,'Twin-Boar')
                        if fuel>=6400*numEng
                            fuel = fuel-6400*numEng;
                        else
                            continue
                        end
                        
                    end
                    
                    obj = RocketStage;
                    obj.mass = Mt+numEng*eng.mass;
                    obj.engines = eng;
                    obj.numEngines = numEng;
                    
                    if densityIdx == 4
                        obj.cost = numEng*eng.cost;
                        fuel = eng.builtInFuel;
                    else
                        %disp(eng)
                        cpp = tankCost./totfuel;
                        cppIDX = knnsearch(totfuel,fuel);
                        CostPerFuel = cpp(cppIDX);
                        %CostPerFuel = interp1(totfuel,tankCost./totfuel,1000, 'nearest','extrap');
                        obj.cost = numEng*eng.cost + CostPerFuel*fuel;
                    end
                    obj.type = 'Linear';
                    obj.dv = stageInfo.dv;
                    obj.fuel = fuel;
                    
                    switch eng.fuelType
                        
                        case 'LF'
                            obj.lFuel = fuel;
                            obj.ox = 0;
                            obj.xenon = 0;
                        case 'LFOX'
                            obj.lFuel = fuel*.45;
                            obj.ox = fuel*.55;
                            obj.xenon = 0;
                        case 'Xenon'
                            obj.xenon = fuel;
                            obj.lFuel = 0;
                            obj.ox = 0;
                        case 'SRB'
                            obj.xenon = 0;
                            obj.lFuel = 0;
                            obj.ox = 0;
                    end
                    newaddition = [stages;obj];
                    outStages = [outStages,newaddition];
                    %if engine meets TWR then dont add more engines
                    break
                end
            end
        end
        
        function outStages = aspStage(stageInfo,stages,engines,tanks)
            %aspStage generates an asparagus staged rocket
            
            g = 9.8;
            %                  LF, LFOX, Xenon, SRB
            FuelDensities = [.005, .005, .001, .0075];
            Ps = .11252; %percent of fueltank that is structural
            %get Mass0
            outStages = [];
            mass0 = 0;
            for i = 1:length(stages)
                mass0 = mass0 + stages(i).mass;
            end
            
            
            %loop through each engine
            switch stageInfo.cond
                
                case 'asl'
                    cond = 1;
                case 'v'
                    cond  = 2;
            end
            
            if stageInfo.e == 0
                %dont loop through engines, because we will just be using
                %the engine from the stage above
                %check that the stage above has an engine
                if stageInfo.center==1
                    error('Tried to make a stage with a center core with no engines, does not make sense')
                end
                if strcmpi(stages(end).engines,'N/A')
                    error('Tried to run an asparagus stage with no engines,\n and the stage above also has no engines');
                end
                
                eng0 = stages(end).engines;
                
                switch eng0.fuelType
                    case 'LF'
                        densityIdx = 1;
                        totfuel = tanks.LF(:,end);
                        tankCost = tanks.LF(:,1);
                    case 'LFOX'
                        densityIdx = 2;
                        totfuel = tanks.RF(:,end) + tanks.RF(:,end-1);
                        tankCost = tanks.RF(:,1);
                    case 'Xenon'
                        densityIdx = 3;
                        totfuel = tanks.Xenon(:,end);
                        tankCost = tanks.Xenon(:,1);
                    case 'SRB'
                        %asparigus staging cannot use SRB
                        outStages = [];
                        return 
                end
                
                numEng0 = stages(end).numEngines;
                eng0.mass = eng0.mass*numEng0;
                T = [eng0.Tasl,eng0.Tv];
                ISP = [eng0.ISPasl,eng0.ISPv];
                T = T.*double(numEng0);
                
                %if there is not enough thrust in the stage above, exit
                if T(cond)<stageInfo.twr*mass0*g
                    outStages = [];
                    return 
                end
                numTankArray = 2*(1:stageInfo.s);
                engArray = repelem(eng0,stageInfo.s);
                
                getAspStageSize = @(x) (stageInfo.dv - RocketStage.aspDV(x,engArray,numTankArray,mass0,cond)).^2;
                %temp = RocketStage.aspDV(1,engArray,numTankArray,mass0,cond)
                %getAspStageSize(1);
                [Mt,value] = fmincon(getAspStageSize,1,[],[],[],[],0);
                if T(cond)<stageInfo.twr*(mass0+Mt*stageInfo.s*2)*g
                    outStages = [];
                    return 
                end
                Mf = (1-Ps)*Mt;
                fuel = Mf/FuelDensities(densityIdx);
                obj = RocketStage;
                obj.mass = Mt*stageInfo.s*2;
                obj.numEngines = 0;
                cpp = tankCost./totfuel;
                cppIDX = knnsearch(totfuel,fuel);
                CostPerFuel = cpp(cppIDX);
                %CostPerFuel = interp1(totfuel,tankCost./totfuel,1000, 'nearest','extrap');
                obj.cost = CostPerFuel*fuel;
                obj.type = ['Asp no Engines S:',num2str(stageInfo.s)];
                obj.dv = stageInfo.dv;
                obj.fuel = fuel;

                switch eng0.fuelType

                    case 'LF'
                        obj.lFuel = fuel;
                        obj.ox = 0;
                        obj.xenon = 0;
                    case 'LFOX'
                        obj.lFuel = fuel*.45;
                        obj.ox = fuel*.55;
                        obj.xenon = 0;
                    case 'Xenon'
                        obj.xenon = fuel;
                        obj.lFuel = 0;
                        obj.ox = 0;
                    case 'SRB'
                        obj.xenon = 0;
                        obj.lFuel = 0;
                        obj.ox = 0;
                end
                newaddition = [stages;obj];
                outStages = [outStages,newaddition];
                %if engine meets TWR then dont add more engines
                return
                
            else
                %loop through engines,
                for eng = engines
                    if stageInfo.center==0
                        eng0 = stages(end).engines;
                        %cost and mass are already reflected in previous
                        %stage
                        eng0.cost = 0;
                        eng0.mass = 0;
                        if isempty(eng0)
                            error('two stages in a row cannot have no engines');
                        end
                        
                    else
                        eng0 = eng;
                    end
                    %if the two engines are not of the same fuel type, exit
                    if ~strcmp(eng.fuelType,eng0.fuelType)
                        continue
                    end
                    
                    switch eng0.fuelType
                        case 'LF'
                            densityIdx = 1;
                            totfuel = tanks.LF(:,end);
                            tankCost = tanks.LF(:,1);
                        case 'LFOX'
                            densityIdx = 2;
                            totfuel = tanks.RF(:,end) + tanks.RF(:,end-1);
                            tankCost = tanks.RF(:,1);
                        case 'Xenon'
                            densityIdx = 3;
                            totfuel = tanks.Xenon(:,end);
                            tankCost = tanks.Xenon(:,1);
                        case 'SRB'
                            %asparigus staging cannot use SRB
                            continue 
                    end
                    for s = 1:stageInfo.s
                        %loop for possible engine numbers
                        for numEng = 1:stageInfo.e
                            if stageInfo.center==0
                                numEng0 = double(stages(end).numEngines);
                            else
                                numEng0 = numEng;
                            end

                            %next step is to make arrays of fake engines

                            numTankArray = [1,1+(2*(1:s))];
                            engArray0 = repelem(eng0,s);
                            masses0 = [eng0.mass*numEng0,repelem(eng.mass*numEng,s)];
                            costs0 = [eng0.cost*numEng0,repelem(eng.cost*numEng,s)];
                            ISPv0 = [eng0.ISPv,repelem(eng.ISPv,s)];
                            ISPasl0 = [eng0.ISPasl,repelem(eng.ISPasl,s)];
                            Tv0 = [eng0.Tv*numEng0,repelem(eng.Tv*numEng,s)];
                            Tasl0 = [eng0.Tasl*numEng0,repelem(eng.Tasl*numEng,s)];
                            
                            engArray = repelem(eng0,s+1);
                            masses = zeros(1,size(engArray0,2));
                            %costs = zeros(1,size(engArray0,2));
                            ISPv = zeros(1,size(engArray0,2));
                            ISPasl = zeros(1,size(engArray0,2));
                            Tv = zeros(1,size(engArray0,2));
                            Tasl = zeros(1,size(engArray0,2));
                            %engArray = [];

                            %make Fake Name For engine Set
                            if stageInfo.center==0
                                centerStr = "c(previous):" + num2str(numEng0) + "x " + eng0.name;
                            else
                                centerStr = "c:" + num2str(numEng0) + "x " + eng0.name;
                            end
                            sizeStr = " S:" + s + " " + num2str(numEng) + "x " + eng.name;


                            for i = 1:size(engArray,2)
                                masses(i) = sum(masses0(1:i));
                                Tv(i) = sum(Tv0(1:i));
                                Tasl(i) = sum(Tasl0(1:i));
                                %Thrust averaged Isp
                                ISPv(i) = sum(ISPv0(1:i).*Tv0(1:i))/Tv(i);
                                ISPasl(i) = sum(ISPasl0(1:i).*Tasl0(1:i))/Tasl(i);
                                engArray(i) = Engine(centerStr+sizeStr,masses(i),Tasl(i),Tv(i),...
                                    ISPasl(i),ISPv(i),0,0,eng.fuelType,0);
                                %disp(engArray0(i))
                            end
                            totCost = sum(costs0);


                            getAspStageSize = @(x) (stageInfo.dv - RocketStage.aspDV(x,engArray,numTankArray,mass0,cond))^2;
                            [Mt,value] = fmincon(getAspStageSize,1,[],[],[],[],0);
                            %test last core stage and launch stage

                            T = [engArray(1).Tasl,engArray(1).Tv];
                            %test core
                            if stageInfo.center==1 && T(cond)<stageInfo.twr * (mass0 + Mt + engArray(1).mass)*g
                                continue
                            end
                            %test launch stage
                            T = [engArray(end).Tasl,engArray(end).Tv];
                            if T(cond)<stageInfo.twr*(mass0 + Mt*numTankArray(end) + engArray(end).mass)*g
                                continue 
                            end

                            Mf = (1-Ps)*Mt;
                            fuel = Mf/FuelDensities(densityIdx);

                            if strcmp(eng.name,'Twin-Boar')
                                if fuel>=6400*numEng
                                    fuel = fuel-6400*numEng;
                                else
                                    continue
                                end

                            end

                            obj = RocketStage;
                            obj.mass = Mt*numTankArray(end) + engArray(end).mass;
                            obj.numEngines = numEng;
                            obj.engines = engArray(end);
                            cpp = tankCost./totfuel;
                            cppIDX = knnsearch(totfuel,fuel);
                            CostPerFuel = cpp(cppIDX);
                            %CostPerFuel = interp1(totfuel,tankCost./totfuel,1000, 'nearest','extrap');
                            obj.cost = CostPerFuel*fuel + totCost;
                            obj.type = 'Asp '+centerStr+sizeStr;
                            obj.dv = stageInfo.dv;
                            obj.fuel = fuel;
                            
                            if isnan(obj.cost)
                                disp('stop here');
                            end
                            
                            switch eng0.fuelType

                                case 'LF'
                                    obj.lFuel = fuel;
                                    obj.ox = 0;
                                    obj.xenon = 0;
                                case 'LFOX'
                                    obj.lFuel = fuel*.45;
                                    obj.ox = fuel*.55;
                                    obj.xenon = 0;
                                case 'Xenon'
                                    obj.xenon = fuel;
                                    obj.lFuel = 0;
                                    obj.ox = 0;
                                case 'SRB'
                                    obj.xenon = 0;
                                    obj.lFuel = 0;
                                    obj.ox = 0;
                            end

                            newaddition = [stages;obj];
                            outStages = [outStages,newaddition];
                            %if engine meets TWR then dont add more engines
                            break
                        end
                    end
                end 
            end  
        end
        
        function dv = aspDV(Mt,engArray,numTankArray,M0,cond)
            %this code pretends that there is one engine that acts as an
            %average of all of the engines combined
            g = 9.8;
            Ps = .11252;
            dv = 0;
            
            %masses = M0;
            for i = 1:length(engArray)
                eng = engArray(i);
                ISP = [eng.ISPasl,eng.ISPv];
                numT = numTankArray(i);
                if i == 1
                    numTEnd = 0;
                else
                    numTEnd = numTankArray(i-1);
                end
                dT = numT-numTEnd;
                mo = M0+Mt*numT+eng.mass;
                mf = M0+eng.mass+Mt*numTEnd+ Mt*dT*Ps;
                %masses = [masses,mf,mo];
                thisDv = ISP(cond)*g*log(mo/mf);
                dv = dv+thisDv;
            end
        end
        
    end
end

