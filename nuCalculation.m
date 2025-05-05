function nu = nuCalculation(sod, xi, Ge, tau0) % вычисление вязкости, где  s0d - sigma0_dimensionless (безразмерный)
    nu = Ge * tau0 * exp(deltaF(sod, xi));
end


