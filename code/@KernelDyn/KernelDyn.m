classdef KernelDyn < KernelBase
% KERNELDYN Kernel-based dynamic programming algorithm for computing an
% optimal control sequence backwards in time. 

    methods
        function obj = KernelDyn(varargin)
        % KERNELDYN Create an instance of the algorithm.

            % Initialization code.
            
            obj = obj@KernelBase(varargin{:});

        end
    end

end
