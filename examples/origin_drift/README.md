#### Example: OriginDrift Analysis
This example demonstrates how to use the generic OriginDrift struct from the moma crate.

It showcases the power of the generic implementation by creating two separate analyzers—one for the PrimeGap strategy and another for the CompositeMass strategy. It feeds both analyzers the same sequence of primes and then compares their final "drift magnitude" scores to determine which strategy produces a more volatile or chaotic sequence of MOMA signatures.

### Concepts Demonstrated
Instantiating OriginDrift with different OriginStrategy types.

Processing a sequence of primes and recording their MOMA signatures.

Calculating the final drift_magnitude to measure volatility.

Comparing the stability of different MOMA strategies.

