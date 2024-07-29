# Handmade SymForce

Port of a chunk of the [SymForce](https://github.com/symforce-org/symforce) 
optimizer to C. Implements just enough to optimize 
[Bundle Adjustment in the Large](https://grail.cs.washington.edu/projects/bal/)
problems.

```
$ /usr/bin/time -l zig-out/bin/balDemo test/bal/data/problem-21-11315-pre.txt
Allocating 1073741824 bytes
BAL optimizer [iter    0] lambda: 1.000000e+00, error prev/new: 4.413239e+06/3.606385e+04, rel reduction: 9.918283e-01
BAL optimizer [iter    1] lambda: 2.500000e-01, error prev/new: 3.606385e+04/3.180932e+04, rel reduction: 1.179722e-01
BAL optimizer [iter    2] lambda: 6.250000e-02, error prev/new: 3.180932e+04/3.108199e+04, rel reduction: 2.286507e-02
BAL optimizer [iter    3] lambda: 1.562500e-02, error prev/new: 3.108199e+04/3.080628e+04, rel reduction: 8.870414e-03
BAL optimizer [iter    4] lambda: 3.906250e-03, error prev/new: 3.080628e+04/3.039174e+04, rel reduction: 1.345638e-02
BAL optimizer [iter    5] lambda: 9.765625e-04, error prev/new: 3.039174e+04/3.037868e+04, rel reduction: 4.297004e-04
BAL optimizer [iter    6] lambda: 2.441406e-04, error prev/new: 3.037868e+04/3.037864e+04, rel reduction: 1.555710e-06
BAL optimizer [iter    7] lambda: 6.103516e-05, error prev/new: 3.037864e+04/3.037864e+04, rel reduction: 9.090088e-09
Success!
nalloc = 0
max_nalloc = 63461280
        0.39 real         0.34 user         0.01 sys
           103612416  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
                6421  page reclaims
                   8  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                   0  signals received
                   2  voluntary context switches
                 686  involuntary context switches
          3947091610  instructions retired
          1090754359  cycles elapsed
           101925760  peak memory footprint
```

## Building from Source

Tested with [Zig](https://ziglang.org/) `0.13.0`

```
zig build
```

Or if you want to go fast,

```
zig build -Doptimize=ReleaseFast
```

