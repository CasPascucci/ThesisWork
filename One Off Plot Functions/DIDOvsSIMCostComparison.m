clear all; close all; clc; format long;

FP2DGCostData = [14.0970, 14.9209, 15.6700, 16.3305, 16.9388, 17.5161, 18.0725, 18.6138, 19.1439, 19.6650, 20.1791, 20.6873, 21.1905, 21.6896, 22.1850, 22.6773, 23.1668, 23.6539, 24.1386, 24.6214, 25.1024]';
DIDOCostData = [13.9769, 14.8818, 15.6382, 16.2894, 16.8909, 17.4630, 18.0155, 18.5540, 19.0817, 19.6010, 20.1136, 20.6206, 21.1228, 21.6211, 22.1157, 22.6073, 23.0962, 23.5827, 24.0670, 24.5492, 25.0298]';

betaVals = [1, .95, .9, .85, .8, .75, .7, .65, .6, .55, .5, .45, .4, .35, .3, .25, .2, .15, .1, .05, 0]';

figure();
plot(betaVals, FP2DGCostData,'r-o');
hold on;
plot(betaVals, DIDOCostData,'b-*');
T = table(betaVals, FP2DGCostData, DIDOCostData, FP2DGCostData - DIDOCostData, ...
      'VariableNames', {'Beta', 'FP2DG', 'DIDO', 'Difference'})