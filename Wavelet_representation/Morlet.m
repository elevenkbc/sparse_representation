function output = Morlet(input)
Sigma = 6;
output = (1/sqrt(1 + exp(-Sigma^2) - 2*exp(-(3/4)*Sigma^2)))*pi^(-1/4)*(exp(-(1/2)*input.^2).*(exp(sqrt(-1)*Sigma*input) - exp(-(1/2)*Sigma^2)));
end