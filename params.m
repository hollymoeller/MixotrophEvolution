function p = params

p(1)  = 0.0;       % a; background mortality (per day)
p(2)  = 1.3;       % b1; bacteria growth rate (per day)
p(3)  = 2.7;       % b2; bacteria growth thermal sensitivity
p(4)  = 50;        % u1; mixotroph growth rate (per day)
p(5)  = 1;         % u2; mixotroph growth sensitivity
p(6)  = 50;        % p1; photosynthetic rate (pg C per cell per day)
p(7)  = 1.88;      % p2; photosynthesis sensitivity
p(8)  = 10;        % g1; grazing rate (bacteria per mixotroph cell per day)
p(9)  = 2.7;       % g2; grazing sensitivity
p(10) = 0.01;      % ra; acclimation respiratory cost (per day)
p(11) = 0.5;       % r0; basal respiration rate (per day)
p(12) = 2.7;       % Q_10,R; respiration thermal sensitivity
p(13) = 0.02;      % lambda; thermal evolution rate (per day)
p(14) = 0.7;       % gamma; N recycling efficiency

p(15) = 10;        % Qminxc (pg per cell)
p(16) = 1.25;      % Qminxn (pg per cell)
p(17) = 0.4;       % bc; bacteria carbon content (pg per cell)
p(18) = 0.1;       % bn; bacteria nitrogen content (pg per cell)

p(19) = 1;         % ki; light half-saturation constant (uE per square meter per day)
p(20) = 1e7;       % kb; bacteria half-saturation constant (cells per mL)
p(21) = 2.0;       % c; plasticity rate (per day)
p(22) = 24;        % T0; reference temperature (deg C)

end