function output = Harr(input)
less_than_half = (0<=input)&(input<0.5);
greater_than_half = (0.5<=input)&(input<1);
output = less_than_half*(1) + greater_than_half*(-1);
end