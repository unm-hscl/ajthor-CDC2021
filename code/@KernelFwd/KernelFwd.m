classdef KernelFwd < KernelBase
% KERNELFWD Kernel based stochastic optimal control algorithm for computing
% an optimal control sequence greedily forward in time. 

    methods
        function obj = KernelFwd(varargin)
        % KERNELFWD Create an instance of the algorithm.

            % Initialization code.
            
            obj = obj@KernelBase(varargin{:});

        end
    end

end
