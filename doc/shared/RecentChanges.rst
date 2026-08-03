.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

The ARKODE Butcher table routine :c:func:`ARKodeButcherTable_CheckOrder` now
generates the Runge-Kutta order conditions programmatically from the rooted
trees associated with the elementary differentials of the ODE right-hand side
(the trees per order follow `OEIS sequence A000081
<https://oeis.org/A000081>`_), raising the maximum analytically-checked order
from 6 to 9 when using double or extended precision. Failed conditions are now
reported using the elementary-weight expression of the corresponding rooted
tree.

**Bug Fixes**

Fixed a missing order-6 condition, ``b'*((A*(c.*c)).*(A*c)) = 1/36``, in
:c:func:`ARKodeButcherTable_CheckOrder`, which could previously over-estimate
the order of a method.

Fixed a bug in :c:func:`ARKodeButcherTable_CheckARKOrder` where the embedding
coefficients of the first table were used in place of those of the second
table when checking the embedding order of an ARK pair.

**Deprecation Notices**
