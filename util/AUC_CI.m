function [CI,SE] = AUC_CI(n_D,n_I,Area)

% By Lukas Krasula
% Inspired by
% *********************  CIAUC  ****************************
%   (c) John W Pickering, Novemeber 2009
%     Christchurch Kidney Research Group
%     University of Otago Christchurch
%     New Zealand
%
%   Last update:  17 July 2012
%
%	Redistribution and use in source and binary forms, with or without 
%   modification, are permitted provided that the following conditions are met:
%
%   * Redistributions of source code must retain the above copyright 
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright 
%     notice, this list of conditions and the following disclaimer in 
%     the documentation and/or other materials provided with the distribution
%
% Attribution to John Pickering.  
% *************************************************************************
% n_D - number of different pairs
% n_I - number of indifferent pairs
% Area - Area under ROC curve

Q1=Area/(2-Area);
Q2=2*Area*Area/(1+Area);

SE=sqrt((Area*(1-Area)+(n_D-1)*(Q1-Area*Area)+(n_I-1)*(Q2-Area*Area))/(n_I*n_D));

CI = 1.96 * SE;