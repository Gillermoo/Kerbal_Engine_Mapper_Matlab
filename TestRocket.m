
clc
clear
close all

%Specify the stages and Allowable Engines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reqs = {'t:pl pl:82.29';
        't:asp e:5 s:3 twr:1.75 cond:asl dv:2800 center:1';
        %'t:asp e:1 s:3 twr:.05 cond:v dv:3000 center:1';
        %'t:asp e:6 s:3 twr:.2 cond:v dv:2700 center:1';
        %'t:lin e:1 twr:.2 cond:v dv:1040';
        %'t:lin e:1 twr:.2 cond:v dv:1660';
        %'t:asp e:6 s:3 twr:1.5 cond:asl dv:2500 center:1';
        %'t:lin e:6 twr:.3 cond:v dv:2500';
        %'t:lin e:6 twr:.3 cond:v dv:2500';
        %'t:lin e:6 twr:.3 cond:v dv:2500';
        %'t:lin e:6 twr:.3 cond:v dv:2500';
        %'t:lin e:6 twr:.75 cond:v dv:2500';
        %'t:lin e:4 twr:.75 cond:v dv:2500';
        %'t:asp e:4 s:4 twr:2 cond:asl dv:2500 center:1';
        %'t:lin e:4 twr:.2 cond:v dv:2300';
        %'t:lin e:4 twr:.6 cond:v dv:1200';
        %'t:asp e:4 s:4 twr:.4 cond:v dv:2400 center:1';
        %'t:asp e:4 s:4 twr:.2 cond:v dv:1700 center:1';
        %'t:pl pl:1.1';
        %'t:lin e:6 twr:.7 cond:v dv:1300';
        %'t:lin e:6 twr:.7 cond:v dv:2600';
        %'t:asp e:4 s:4 twr:.7 cond:v dv:2900 center:1';
        %'t:lin e:4 twr:.7 cond:asl dv:2900'
        %'t:asp e:4 s:4 twr:1.25 cond:asl dv:2500 center:1'
        %'t:lin e:6 twr:1.25 cond:asl dv:2500';
        };
%{   
't:pl:5';
't:lin e:4 twr:.2 cond:v dv:1000';
't:asp e:4 s:4 twr:1.25 cond:asl dv:3600 center:1'
%}
    
Engines = {
'Spider'        ,1;
'Twitch'        ,1;
'Thud'          ,1;
'Ant'           ,1;
'Spark'         ,1;
'Terrier'       ,1;
'Reliant'       ,1;
'Swivel'        ,1;
'Vector'        ,1;
'Dart'          ,1;
'Nerv'          ,0;
'Poodle'        ,1;
'Skipper'       ,1;
'Mainsail'      ,1;
'Twin-Boar'     ,1;
'Rhino'         ,1;
'Mammoth'       ,0;
'R.A.P.I.E.R.'  ,1;
'Dawn'          ,1;
'Mastodon'      ,0;
'Cheetah'       ,1;
'Bobcat'        ,1;
'Skiff'         ,1;
'Wolfhound'     ,0;
'Kodiak'        ,1;
'Cub'           ,1;
'Flea'          ,1;
'Hammer'        ,1;
'Thumper'       ,1;
'Kickback'      ,1;
'Sepratron I'   ,1
'Shrimp'        ,1;
'Mite'          ,1;
'Thoroughbred'  ,1;
'Clydesdale'    ,1;
'Pullox'        ,1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('LFTanks.mat');
load('RFTanks.mat');
load('XenonTanks.mat');
tanks.LF = LFTanks;
tanks.RF = RFTanks;
tanks.Xenon = XenonTanks;

allowEngines = Engine.setupEngines(Engines);
r1 = Rocket(reqs,allowEngines,tanks);
%sound(sin(linspace(0,1000,5000)), 8192)

%s1 = RocketStage();

maxdv = 2800;
span = 25;
dv1 = linspace(0,maxdv,span);
dv2 = maxdv-dv1;

%master = 't:asp center:1 e:5 s:3 twr:.3 cond:v dv:';
master = 't:lin e:5 twr:1.75 cond:asl dv:';
type = 'cost';

outlinegeoquadVAL = [];
%outrockets = [];

for i = 1:span
    reqs2 = reqs;
    if(dv1(i)~=0)
        reqs2 = [reqs2;cat(2,master,num2str(dv1(i)))];
    end
    if(dv2(i)~=0)
        reqs2 = [reqs2;cat(2,master,num2str(dv2(i)))];
    end
    try
        r1 = Rocket(reqs2,allowEngines,tanks);
        if(strcmp(type, 'cost'))
            vals = {r1.cost};
            [val,idx] = min([vals{:}]);
            r0 = r1(idx);
            outval(i) = r0.cost;
            outrockets(i) = r0;
        else
            vals = {r1.mass};
            [val,idx] = min([vals{:}]);
            r0 = r1(idx);
            outval(i) = r0.mass;
            outrockets(i) = r0;
        end
    catch
        outval(i) = nan;
        outrockets(i) = Rocket.empty;
    end
    
    
end

plot(dv1,outval);
[val,idx] = min(outval);
r1 = outrockets(idx);

%a = repelem(s1,10)

%r1 = Rocket([s1,s1,s1]);
%r1.toString();










