%Author: William Craver
%Date:   11/7/2020
%Desription: this file sets up which engines you want to plot and plots the
%optimal cost and mass plots for the given engines.  All engines are listed
%here so you can uncheck the ones you might not have becasue of DLC or not
%unlocked yet.

clc
clear 
close all

m = logspace(-1,3,800);
dv = linspace(100,6200,800);
e = 1;

[DV,M] = meshgrid(dv,m);

info.t = 'lin';
info.e = 1;
info.twr = 2.5;
info.cond = 'asl';
info.dv = DV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
'Nerv'          ,1;
'Poodle'        ,1;
'Skipper'       ,1;
'Mainsail'      ,1;
'Twin-Boar'     ,1;
'Rhino'         ,1;
'Mammoth'       ,1;
'R.A.P.I.E.R.'  ,1;
'Dawn'          ,1;
'Mastodon'      ,1;
'Cheetah'       ,1;
'Bobcat'        ,1;
'Skiff'         ,1;
'Wolfhound'     ,1;
'Kodiak'        ,1;
'Cub'           ,1;
'Flea'          ,1;
'Hammer'        ,1;
'Thumper'       ,1;
'Kickback'      ,1;
'Sepratron I'   ,1;
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
conds = {'v'};
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'Vac.gif'; 

 % Capture the plot as an image 

h = RocketStage.plotLinearStage(info,M,allowEngines,tanks);
    
     







