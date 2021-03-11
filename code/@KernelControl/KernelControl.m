classdef KernelControl < KernelBase
% KERNELCONTROL Kernel based stochastic optimal control algorithm for
% computing the one step optimal control actions.

    methods
        function obj = KernelControl(varargin)
        % KERNELCONTROL Create an instance of the algorithm.

            % Initialization code.
            
            obj = obj@KernelBase(varargin{:});

        end
    end

end
