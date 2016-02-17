Tilings of non-face-aligned cubes around the origin
===================================================

This is some (somewhat ugly but fast) C code to try to solve the problem.

It could be better documented if requested.

It's a work in progress.

layout
------
TODO: brute failed because there are just too many possible arrangements, there are far fewer layouts.  Rather than computing all arrangements and classifying the histograms, let's take each arrangement and see if it's possible.

This is making good progress.  I've put some terrible hacks in it; there are lots of bit fields floating around that make for very unintuitive code, but they a major speedup.  There's also some terrible hackery with includign buildlayout.c multiple times with different #defines; this avoids some conditionals that actually result in a significant overall speedup.  Sadly, with code this tight, those conditionals actually matter (my unconfirmed suspicion is that adding them in results in extra register spillage, or maybe they're poorly predicted).  These hacks have given a 33% speedup over already very fast code.

### Pruning by matchings
One big win is finding the size of the maximal matching remaining, and then ensuring that it's roughly the right size.  More specifically:

Consider the number of cubes left that must be matched, either as part of a 1-dimensional pairing or a $d$-dimensional $2^d$ set.  To pair in the former case, this requires a pair of adjacent cubes; in the latter, $2^{d-1}$ pairs of adjacent cubes.  In either case, we can compute how many cubes must be paired to complete the histogram.  (Or equivalently, we know how many 0-dimensional cubes must remain unpaired at the end; the rest of the cubes must be paired.)

Since the implied graph over the hypercube is bipartite (give each vertex of the hypercube a coordinate from {0, 1}^n, and split into odd and even parities), we can compute the maximum matching over the remaining unmatched vertices using standard techniques (ideally, tweaked slightly to improve performance).  If that maximal matching is too small to cover the remaining pairings in the histogram, we needn't bother searching it - there are too many isolated cubes to allow us to pair enough up.  The inverse doesn't work - there may be numerous matchings of the right size, but they may be face-aligned, so we do need to search in that case.  But this can let us prune huge numbers of cases with few unmatched cubes, reducing runtime from decades to minutes.

Moreover, this demonstrates that any histogram with more than half its cubes unpaired can't be made.

It would also be nice to be able to find a minimum maximal matching - a matching that can't be extended, but is as small as possible.  This is because a non-face-aligned layout can't have any unpaired adjacent 0-dimensional cubes; they would be trivially face-aligned.  So if the minimal maximal matching is larger than the number of pairings remaining, then no matter how those pairings are chosen, there will remain an unpaired adjacent pair.  (E.g., suppose we have five one-dimensional pairings left to place, but the minimum maximal matching is of size six.  Then no matter how we place those five pairings, they can't f orm a maximal matching, so that placement can be extended with a sixth pairing, which represents two adjacent unpaired (and hence, face-aligned) cubes.)

Unfortunately, finding a minimal maximal matchings is NP-hard.  But simply taking half the size of the maximum matching provides a lower bound - consider overlaying a matching with $k/2 - 1$ or fewer edges over a maximum matching with $k$ edges.  Each of the $\le k/2 - 2$ edges of the matching can cover at most two vertices used by the maximum matching, for a total of at most $k - 2$ vertices covered.  But we can now pick the two vertices in the maximal matching uncovered by the smaller matching, and use those as the ends of an augmenting path to extend the matching, so it can't be maximal.  (This is a little fuzzy, but Wikipedia confirms that the maximal matching is a 2-approximation of the minimum maximal matching.)

Or, for a simpler (and less handwavy) proof, consider the $k$ edges of the maximum matching directly.  Each of the $< k/2$ edges of a hypothetical smaller matching covers two vertices; these two vertices are members of at most two edges of the maximum matching.  So there are less than $2(k/2) = k$ edges with at least one vertex covered by the smaller matching, so there is not only an augmenting path, but a simple edge of the maximum matching which is unmatched.  Clearly, this smaller matching can't be maximal.

So we can simply compute the maximal matching, and check that the number of remaining pairings to be done is no more than that number, and at least half of it.  This is a *huge* win, reducing the number of unknown histograms from 249 to 98 in under a week, without a significant amount of optimization effort.

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
