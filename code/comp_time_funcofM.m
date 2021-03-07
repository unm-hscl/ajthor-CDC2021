classdef comp_time_funcofM < matlab.perftest.TestCase
% COMP_TIME Performance testing class to measure computation time of the
% algorithm.
    
    properties
        
        % Number of test points.
        R = 10;
        
        % Variables to hold samples.
        xs
        us
        ys
        
        c
        
        % Algorithm instantiation.
        alg
        
    end
    
    properties (MethodSetupParameter)
        
        % Range of the number of test points to test. 
        NumberOfSamplePoints = num2cell(100:500:8100);
        
    end
    
    methods (TestMethodSetup)
        
        function setupAlgorithm(testCase, NumberOfSamplePoints)
            
            % Sampling time.
            Ts = 0.25;
            
            % System matrices.
            A = [1 Ts; 0 1];
            B = [(Ts^2)/2; Ts];
            
            % Samples.
            testCase.xs = [
                -1 + 2*rand(1, NumberOfSamplePoints);
                -1 + 2*rand(1, NumberOfSamplePoints);
                ];

            testCase.us = -1.1 + 2.2*rand(1, size(testCase.xs, 2));
            
            W = 0.01*randn(size(testCase.xs));

            testCase.ys = A*testCase.xs + B*testCase.us + W;
            
            % Evaluate the cost function at the sample points.
            testCase.c = vecnorm(testCase.ys);
            
            % Instantiate the algorithm.
            testCase.alg = KernelControl('Lambda', 1/NumberOfSamplePoints^2, ...
                                            'Sigma', 0.25);
            
        end
        
    end
    
    methods (Test)
        
        function testAlgorithmVaryingTestPoints(testCase)

            % Generate test points.
            rng(0);
            xt = rand(2, 10);
            
            % Admissible control inputs.
            ur = linspace(-1, 1, 100);
            
            testCase.startMeasuring();
            
            testCase.alg.compute(testCase.xs, testCase.us, ...
                                    xt, ur, ...
                                    testCase.c);
                                
            testCase.stopMeasuring();
            
        end
        
    end
end