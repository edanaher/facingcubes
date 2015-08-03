Tilings of non-face-aligned cubes around the origin
===================================================

This is some (somewhat ugly but fast) C code to try to solve the problem.

It could be better documented if requested.

It's a work in progress.

layout
------
TODO: brute failed because there are just too many possible arrangements, there are far fewer layouts.  Rather than computing all arrangements and classifying the histograms, let's take each arrangement and see if it's possible.

This is making good progress.  I've put some terrible hacks in it; there are lots of bit fields floating around that make for very unintuitive code, but they a major speedup.  There's also some terrible hackery with includign buildlayout.c multiple times with different #defines; this avoids some conditionals that actually result in a significant overall speedup.  Sadly, with code this tight, those conditionals actually matter (my unconfirmed suspicion is that adding them in results in extra register spillage, or maybe they're poorly predicted).  These hacks have given a 33% speedup over already very fast code.


brute
-----
First attempt - start with $2^n$ individual cubes at $\{\pm 1}^n$, then iteratively "merge" adjacent pairs, replacing them with a single cube with coordinate 0 along the shared axis.  This is in two steps: first, find a matching of the first layer of cubes that leaves no cubes at the first layer sharing a face; second, merge these matched cubes (resulting in cubes with two zero coordinates), and repeat until finished.  This layer-by-layer merging makes it easier to avoid face-aligned cubes; these are simply cubes that are unmatched when their layer is completed.  (And the two-phase process takes advantage of the fact that the initial layer is extremely dense, having all $2^n$ cubes, but further cubes are only merged when they share the same set of zero dimensions, so these sets are very sparse.) 

There's a good amount of C cleverness - the coordinates are stored as packed bits in an integer, with 1 and 0 representing 1 and -1, with merged zero coordinates stored out-of-band.  This allows from some nice bit operations; e.g., the "merged" cube is simply the binary and of its two inputs, and checking if two cubes are adjacent is simply `!((a^b) & ((a^b)-1)`.

Between these tricks and a number of ["don't do anything stupid"/"don't do anything twice"](http://www.dcc.fc.up.pt/~pribeiro/estagio2008/usaco/4_1_Optimization.htm) optimizations, finding all arrangements of 5-dimensional cubes takes just over a second on a ~4 year old desktop.  Unfortunately, 6 dimensions is infeasible; for 5 dimensions, there were about a million final solutions checked.  For 6 dimensions, the current program's search needs to check about 35 trillion arrangements only for the first 16 matched pairs of the top layer, and each extra pair is adding roughly a factor of 3 at this point.  A 9-hour overnight run fully checked 2300 arrangements of 21 pairs; this suggests that the total run time would be $35e12 \times 3^5 \times 9 \approx 76e15$ hours, or over 8 trillion years.

Granted, this approach should parallelize nearly perfectly, and there's probably still a bit more optimization to be done, so it's probably possible to do this within the lifetime ofthe universe.  But still, this approach seems unlikely to produce a solution within a human lifetime.

Results
-------
Each of the following files is also available in sorted form, which sorts the histograms in a more convenient manner:

- results/brute.\[n\].out - For 3, 4, and 5 dimensions, these files are the final output.  The first column is the number of times each given histogram was found.  The $k$th value histogram itself is a count of how many cubes have $k$ zeros in their coordinates; equivalently, the first column is the number of unmerged cubes, the second is the number of pairs of cubes merged once, and so on, with the last column being the number of fully-merged cubes (obviously always zero, except for the trivial single cube covering).

- results/brute.6.wip - For 6 dimensions, a complete run is infeasible.  But I've been running it for over three weeks now to see how it does; this is its output so far.

- results/random.6.wip - This run did a random selection of the search tree.  After 5 runs (a couple days, IIRC), it found a good number of solutions, many of which are disjoint from what the full search has found so far.

I'm wary of estimating the total number of solutions from the combination of brute.6 and random.6, since I'd expect there to be dramatic nonuniformity in both of them, making extrapolating seem difficult (if not impossible) without a deeper understanding of the structure of the problem.

- counts.6.d.17: In an attempt to quantify how much caching helped, these are the number of arrangements checked at up to 17 pairings, caching up to depth d.  Compared to d=0 (no caching), there is a ~40x speedup with diminishing returns for caching at deeper levels.  This also demonstrates the dramatic blowup that suggests the impracticality of this approach to the problem.
