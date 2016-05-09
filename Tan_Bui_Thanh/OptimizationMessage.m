function OptimizationMessage(exitflag,OutputFile)
% OPTIMIZATION MESSAGE
% Output message for the optimization solver
% Copyright (c) 2007 by Tan Bui
% All Rights Reserved. Massachusetts Institute of Technology

Message={'Cost is below the tolerance','Step is less then the tolerance',...
    'TR size is too small','Exceeds MaxIter',...
    'Scaled gradient is less than the tolerance','STIRNCG fails',...
        'Solution blows up', 'shape parameters are not feasible'};

if nargin < 2
  fprintf('%s\n',char(Message(exitflag)));
else
  fprintf(OutputFile,'%s\n',char(Message(exitflag)));
  
end
