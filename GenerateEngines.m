%Author: William Craver
%Date:   11/7/2020
%Desription: this script creates the .mat files that are loaded to run the
%rocket sims,

clc
clear
close all

[num, txt, raw] = xlsread('LFTanks.csv');

CostPerFuel = num(:,1)./num(:,7);
p = paretoPoints([num(:,7),CostPerFuel],[-1,-1]);

LFTanks = num(p(:,end),:);
[~,idx] = sort(LFTanks(:,end));
LFTanks = LFTanks(idx,:);
LFTanks(3,:) = [];
save('LFTanks.mat','LFTanks');

[num, txt, raw] = xlsread('RFTanks.csv');

CostPerFuel = num(:,1)./num(:,7);
p = paretoPoints([num(:,7),CostPerFuel],[-1,-1]);

RFTanks = num(p(:,end),:);
[~,idx] = sort(RFTanks(:,end));
RFTanks = RFTanks(idx,:);
save('RFTanks.mat','RFTanks');

[num, txt, raw] = xlsread('XenonTanks.csv');

CostPerFuel = num(:,1)./num(:,7);
p = paretoPoints([num(:,7),CostPerFuel],[-1,-1]);

XenonTanks = num(p(:,end),:);
[~,idx] = sort(XenonTanks(:,end));
XenonTanks = XenonTanks(idx,:);
save('XenonTanks.mat','XenonTanks');




