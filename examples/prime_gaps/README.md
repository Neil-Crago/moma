# Example: Prime Gap Analyzer

This example demonstrates how to use the analytical tools within the `moma` crate to study the statistical properties of the gaps between consecutive prime numbers.

It showcases the `PrimeGapField` struct to calculate barycentric offsets, modular class distributions, and entropy scores, providing insights into the structure and randomness of the prime sequence.

-----

## Concepts Demonstrated

  * Creating a `PrimeGapField` from a sequence of primes.
  * Calculating the **barycentric offset** of each gap (its deviation from the local average).
  * Analyzing the distribution of prime gaps modulo `n` (in this case, `mod 6`).
  * Calculating the **Shannon entropy** for each modular class to measure the randomness of the distribution.
  * Filtering gaps to find statistical **outliers**.
  * Using the collected prime data to project potential **Goldbach pairs**.

-----

## How to Run

From the root of the `moma` crate repository, run the following command:

```sh
cargo run --example prime-gaps
```

-----

## Expected Output

You should see output similar to the following, which summarises the analysis of primes up to 100:

```
--- Prime Gap Field Analysis --- 🌌

📊 Entropy Scores for Prime Gaps mod 6:
   Class  0: 0.2575
   Class  2: 0.5288
   Class  4: 0.5211

🌠 Outlier Gaps (Barycentric Offset > 3.0):
   Gap [ 89,  97] | Size:  8 | Offset: +3.67

✨ Goldbach Projections for 96:
   96 = 5+91 = 7+89 = 13+83 = 17+79 = 19+77 = 23+73 = 29+67 = 37+59 = 43+53
```