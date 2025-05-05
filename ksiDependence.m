function ksi = ksiDependence(eps0)
    ksi = (260 + eps0 .* (123200 + (76250-114200 * eps0) .* eps0))./(12400 + eps0 .* (16770 + eps0 .* (-20670 + (-8820 + eps0) .* eps0)));
end

