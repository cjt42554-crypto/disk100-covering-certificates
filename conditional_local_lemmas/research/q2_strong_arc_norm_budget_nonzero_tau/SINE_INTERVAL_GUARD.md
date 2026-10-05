# Endpoint Taylor guard for sine on [1.5,1.65]

The frozen full proof's needed estimate is `sin z>.996` on
`1.5<=z<=1.65`. It may be proved without a precise enclosure of pi.
Sine is concave on this interval, which lies inside `(0,pi)`.
Therefore its minimum is attained at an endpoint.
Use the seventh-order lower alternating Taylor polynomial

    P7(z)=z-z^3/6+z^5/120-z^7/5040 <=sin z.

At z equal to 1.5, the first three terms total `1.00078125`,
and `z^7/5040<.0034`, since
`17.0859375<(.0034)(5040)=17.136`. Hence

    sin(1.5)>1.00078125-.0034=.99738125>.996.

At z equal to 1.65, the first three terms total
`1.0032275859375`, and `z^7/5040<.00661`, since
`33.29565857578125<(.00661)(5040)=33.3144`. Hence

    sin(1.65)>1.0032275859375-.00661
              =.9966175859375>.996.

Concavity proves the strict whole-interval estimate. This replaces
the optional distance-to-pi/2 argument for that guard only. The
full proof's source and all earlier manifests remain unchanged.
All products above are paper exact decimals; mathematical program
execution, imports, checker replay, scan and search: 0.
