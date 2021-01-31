function S=StartingDesign(initial_design,number_startpoints, Data)

% generates initial experimental design in one of four possible ways
%--------------------------------------------------------------------------
%Copyright (c) 2012 by Juliane Mueller
%
% This file is part of the surrogate model module toolbox.
%
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%Tampere University of Technology, Finland
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%input:
%initial_design - string, name of design strategy 
%number_startpoints - integer, number of points desired for the initial experimental design
%Data - structure, contains all problem information
%
%output:
%S - matrix with design points, (m x d), m=number of points, d = dimension
%--------------------------------------------------------------------------
%LHS - matlab's integrated Latin hypercube sampling (maximizes minimum
%distance)
%SLHD - symmetric Latin hypercube design
%SPACEFIL - startegy as used in EGO (see Forrester implementation
%CORNER - uses all or several corner points plus mid point 
%--------------------------------------------------------------------------

if strcmp(initial_design,'LHS')%lhsdesign
    S = lhsdesign(number_startpoints,Data.dim,'criterion','maximin',...
            'iterations',20);
    S=repmat(Data.xup-Data.xlow,number_startpoints,1).*S+repmat(Data.xlow,number_startpoints,1);    
elseif strcmp(initial_design,'SLHD')%%symmetric latin hypercube
    S=SLHD(Data.dim,number_startpoints);
    S=repmat(Data.xup-Data.xlow,size(S,1),1).*S...
        +repmat(Data.xlow,size(S,1),1);    
elseif strcmp(initial_design,'SPACEFIL')%spacefilling design as in EGO
    S=bestlh(number_startpoints,Data.dim,20,10);
    S=repmat(Data.xup-Data.xlow,size(S,1),1).*S...
        +repmat(Data.xlow,size(S,1),1);    
elseif strcmp(initial_design,'CORNER') %corner point strategy
    display('Corner point strategy may fail, still in beta')
    S=cornerpoints(Data, number_startpoints);
end
    
end %function