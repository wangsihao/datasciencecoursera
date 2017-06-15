function approx = get_function(params, t)

x = length(params);

if (x == 3)
    approx = params(1) + params(2)*exp(-t/params(3));
else
    approx = params(1) + params(2)*(exp(-t/params(3)) - exp(-t/params(4)));
end