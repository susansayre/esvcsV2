function [basisCoefficients,basis] = approximateFcn(basis,nodeVals,func,P)

basisCoefficients = broyden(func,zeros(n,1),basis,P);