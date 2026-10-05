# Cycle-coordinate decision

## Archive finding

The older cyclic-cosmos program names three coordinates:

- local oscillator phase theta;
- cosmic scalar coordinate phi;
- global cycle parameter lambda_cycle.

The archive contains working descriptions of theta and phi, but the search for a mathematical definition of lambda_cycle returns only the name in the sector list; no reproducible formula or approved normalization is supplied.

## Production decision

The production engine therefore implements:

- local phase theta from the solved local S/Sdot state;
- cosmic scalar phi from the solved COSMOS scalar;
- H_eff crossing events from the solved unified geometry;
- local proper time from the solved central lapse;
- derived e-fold coordinate from integral H_eff dt.

lambda_cycle remains undefined until the archive supplies an actual equation or independently derived definition.

This is intentional. A synthetic cycle parameter would be a hidden imposed coordinate, not an observed quantity.
