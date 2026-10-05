# Constructive covers and exact finite comparisons

This directory contains selected mathematical data and verification sources. All reported mathematical PASS results were completed before this public preparation and are explicitly inherited. Public preparation copied data, adapted paths, parsed source syntax and bound hashes; it ran no covering verifier, Hall proof, optimizer or search. The two new portable adapters have not been mathematically executed.

The guarded endpoint consists of 100 actual binary64 centers. Its exact finite-decimal CSV represents those same 200 coordinates. The inherited integer/surd verifier proved a covering of the closed unit disk at radius `0.116404851624 = 14550606453/125000000000`. This is a configuration-specific upper bound and is larger than the project's record upper bound; it is not an optimality proof.

The saved ordinary endpoint is within `delta = 1/100000000` of the guarded endpoint center by center. All 100 squared displacement checks used exact dyadic arithmetic. The triangle inequality therefore transfers the covering to radius `0.116404861624`. It also transfers the guarded endpoint's finite comparison bound `0.002362504` to the ordinary endpoint's strict bound `0.002362494` against the same references.

The Hall certificate uses exactly 35,599 saved reference arrays, each with all 100 indexed centers. For every reference, every orthogonal transform fixing the origin and every bijection of the 100 centers has maximum matching distance strictly greater than `0.002362504`. Rotations and reflections are both included. Translation and scaling are outside the statement. The references are a finite saved set; this does not exclude all configurations, all historical arrays, a basin or a global optimum. The guarded endpoint itself is not one of these references.

For rational `t >= 0` and `epsilon > 0`, let `L` count centers of one array with norm at most `t`, and let `H` count centers of the other with norm at most `t + epsilon`. If `L > H`, no matching with maximum distance at most `epsilon` exists: norm preservation and the reverse triangle inequality would inject the first indexed set into the second. Equality is included in both threshold sets. Coincident centers retain separate indices.

The unchanged `common/exact_radial.py` acceptance kernel represents each squared norm as `N/2^q`. For `t = u/D >= 0`, membership is the integer comparison `N D^2 <= u^2 2^q`. Each saved witness contains a rational threshold, both 100-bit masks, their exact counts, and a second buffered pair. A buffer `b` replaces the two thresholds by `t-b` and `t+epsilon+b`, proving the stronger bound `epsilon+2b`. The minimum exact buffer supplies the common finite bound. Numeric ranks and gaps in the saved witness file were generation diagnostics and do not enter the public reader's acceptance.

The exact cover verifier checks all interior circumcenter candidates, all pair-bisector/unit-circle intersections, each nonzero center's antipode, and the four cardinal points. These contain a maximum of nearest-center distance on the closed disk: an interior maximum occurs at a Voronoi vertex, and a boundary maximum occurs at a Voronoi break or a single-center boundary extremum. Candidate inclusion and disk coverage use exact integers or exact quadratic-surd sign comparisons. Floating calculations only select a disk witness and report approximations.

From the repository root, optional future verification commands are:

```sh
python3 -S -B constructive/endpoint/source/verify_cover.py constructive/endpoint/inputs/centers.csv 0.116404851624
python3 -S -B constructive/hall/source/check_all.py
python3 -S -B constructive/transport/source/verify_transport.py
```

The cover verifier has a standard-library fallback if SciPy is unavailable. The new Hall and transport adapters use only the Python standard library. They resolve paths relative to this directory and verify the clean manifest and input hashes. The transport command checks exact displacement; its transferred statements also require the separate source cover and Hall certificates. No command is run automatically.

The two `original_generator.py` files preserve the original mathematical producer sources byte for byte as source records. Their archival launch/admission interfaces and floating selector dependencies are not part of the portable public entrypoints. Use the commands above. `SOURCE_PUBLIC_CORRESPONDENCE.json` identifies unchanged data and kernels, the compact inventory transformation and the new adapter mappings using project-relative names and hashes only.

All 36 coordinate chunks and the saved witness file retain their original mathematical bytes. The clean inventory retains the exact reference entry order, chunk rows and coordinate digests, while omitting historical aliases, personal paths and operational metadata. The 122 unrecovered prior ordered hashes and 122 unrecovered prior unordered hashes are not coordinate-backed by this finite certificate and remain outside its conclusion.
