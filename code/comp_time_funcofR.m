classdef comp_time_funcofR < matlab.perftest.TestCase
% COMP_TIME Performance testing class to measure computation time of the
% algorithm.
    
    properties
        
        % Number of sample points.
        M = 5625;
        
        % Variables to hold samples.
        xs
        us
        ys
        
        c
        
        % Algorithm instantiation.
        alg
        
    end
    
    properties (TestParameter)
        
        % Range of the number of test points to test. 
        NumberOfTestPoints = num2cell(1:5:51);
        
    end
    
    methods (TestMethodSetup)
        
        function setupAlgorithm(testCase)
            
            % Sampling time.
            Ts = 0.25;
            
            % System matrices.
            A = [1 Ts; 0 1];
            B = [(Ts^2)/2; Ts];
            
            % Samples.
            testCase.xs = [
                -1 + 2*rand(1, testCase.M);
                -1 + 2*rand(1, testCase.M);
                ];

            testCase.us = -1.1 + 2.2*rand(1, size(testCase.xs, 2));
            
            W = 0.01*randn(size(testCase.xs));

            testCase.ys = A*testCase.xs + B*testCase.us + W;
            
            % Evaluate the cost function at the sample points.
            testCase.c = vecnorm(testCase.ys);
            
            % Instantiate the algorithm.
            testCase.alg = KernelControl('Lambda', 1/testCase.M^2, ...
                                            'Sigma', 0.25);
            
        end
        
    end
    
    methods (Test)
        
        function testAlgorithmVaryingTestPoints(testCase, NumberOfTestPoints)

            % Generate test points.
            rng(0);
            xt = rand(2, NumberOfTestPoints);
            
            % Admissible control inputs.
            ur = linspace(-1, 1, 100);
            
            testCase.startMeasuring();
            
            testCase.alg.compute(testCase.xs, testCase.us, ...
                                    xt, ur, ...
                                    testCase.c);
                                
            testCase.stopMeasuring();
            
        end
        
        function testOptimalAlgorithmVaryingTestPoints(testCase, NumberOfTestPoints)

            % Generate test points.
            rng(0);
            xt = rand(2, NumberOfTestPoints);

            % Sampling time.
            Ts = 0.25;

            % System matrices.
            A = [1 Ts; 0 1];
            B = [(Ts^2)/2; Ts];
            
            dynamics = @(x, u) A*x + B*u;

            % Compute the results using the optimization based algorithm.
            alg_opt = OptimalControl();
            
            testCase.startMeasuring();
            
            alg_opt.compute(dynamics, xt, [-1, 1]);
            
            testCase.stopMeasuring();

        end
        
    end
end