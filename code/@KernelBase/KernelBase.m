classdef KernelBase < handle
%KERNELBASE Abstract base class for kernel algorithms.

    properties

        % Algorithm parameters.
        Lambda (1, 1) double = 1
        Sigma (1, 1) double = 0.1

    end
    
    methods
        function obj = KernelBase(varargin)
        % KERNELCONTROL Create an instance of the algorithm.

            % Initialization code.
            
            p = inputParser;
            addParameter(p, 'Lambda', 1, ...
                @(arg) validateattributes(arg, {'numeric'}, {'nonnegative'}));
            addParameter(p, 'Sigma', 0.1, ...
                @(arg) validateattributes(arg, {'numeric'}, {'nonnegative'}));
            parse(p, varargin{:});
            
            obj.Lambda = p.Results.Lambda;
            obj.Sigma = p.Results.Sigma;
            
        end
    end
end

