% Function ModelCenterWrap(x,fidelity):
%   Calls a java program (compiles if necessary) which instantiates the
%   ModelCenter API, loads a ModelCenter program that is set in the java
%   program and outputs the result to the command line.
%
%   This wrapper calls the java program, reads the result from the command
%   line and returns it.
%
%   Inputs:
%       x - design vector, a column vector of inputs to run
%       fidelity - an integer fidelity flag sent to ModelCenter
%
%   Outputs:
%       f - the scalar objective function value from ModelCenter
%
function f = ModelCenterWrap2(x,fidelity,varargin)

% Note this if statement keeps the low-fidelity model in Matlab and the
% high-fidelity model in ModelCenter. To have both models in ModelCenter
% simply remove the if statement. The final parameter passed to ModelCenter
% is the fidelity level, it is passed regaurdless of where the low-fidelity
% model is located.
if(fidelity==2)
    % Path to the Java Development Kit:
    path_to_java='c:\Progra~1\Java\jdk1.6.0_21\bin\';

    str=cd;
    % Path to the Main.java file to run your ModelCenter model:
    cd('C:\Program Files\Phoenix Integration\ModelCenter 9.0');

    % Call to compile Java code if necessary:
    [r,c]=system([path_to_java,'javac -classpath "c:\Program Files\Phoenix Integration\ModelCenter 9.0\ModelCenter.jar" Main.java']);

    argString = num2str(x');
    fidString = num2str(fidelity);

    % Call to Actually run your ModelCenter Model:
    [r,c] = system(['java -Djava.library.path="C:\Program Files\Phoenix Integration\ModelCenter 9.0" -classpath "c:\Program Files\Phoenix Integration\ModelCenter 9.0\ModelCenter.jar;c:\Program Files\Phoenix Integration\ModelCenter 9.0" Main ',argString,' ',fidString]);

    % Convert ModelCenter output to double (if not already):
    if(ischar(c))
        f=sscanf(c,'%f');
    else
        f=c;
    end
    % Back to working directory:
    cd(str);
elseif(fidelity==1)
    f=x'*x;
else
    f=0;
end