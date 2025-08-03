# Example: MOMA for Bioinformatics

This example demonstrates a creative application of the MOMA framework to bioinformatics. It explores the idea of mapping abstract mathematical patterns derived from prime numbers to concrete biological events, such as genetic mutations, within a simulated chaotic environment.

## Concepts Demonstrated

* **Interdisciplinary Application**: Bridging number theory (MOMA) and bioinformatics.
* **Event Simulation**: Using a numeric MOMA signature to determine the location of a simulated point mutation in a DNA sequence.
* **Mutation Analysis**: Classifying the effect of each mutation (Silent, Missense, or Nonsense) based on the standard genetic code.
* **Entropy Pulses**: Tracking the Shannon entropy of recent signatures to model environmental chaos. A spike in entropy is treated as a "chaotic event" or "entropy pulse," simulating a stressor that might increase mutation rates.
* **Data Visualization**: Using the `plotters` crate to create a heatmap of the generated signatures, providing a unique look at the structure of the chosen MOMA strategy.

## How to Run

First, ensure you have the `plotters` dependency in your main `Cargo.toml`:

```toml
[dependencies]
plotters = "0.3.5"
```
When you run the simulation it will print out each structure-forming event and any subsequent chaotic eras that are triggered, and will provide a final summary report.