classdef KernelControl < handle
% KERNELCONTROL An algorithm.

    properties (Access = private)

        % Algorithm parameters.
        lambda_ (1, 1) double = 1
        sigma_ (1, 1) double = 0.1

    end

    methods
        function obj = KernelControl(varargin)
        % KERNELCONTROL Create an instance of the algorithm.

            % Initialization code.
            
            p = inputParser;
            addParameter(p, 'Lambda', 1, @(arg) validateattributes(arg, {'numeric'}, {'nonnegative'}));
            addParameter(p, 'Sigma', 0.1, @(arg) validateattributes(arg, {'numeric'}, {'nonnegative'}));
            parse(p, varargin{:});
            
            obj.lambda_ = p.Results.Lambda;
            obj.sigma_ = p.Results.Sigma;

        end
    end

end
