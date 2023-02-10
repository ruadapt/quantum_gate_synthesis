# quantum_gate_synthesis

This project is a C++ implementation of parts of [newsynth](https://hackage.haskell.org/package/newsynth), a Haskell library. Most of the files in the Haskell implementation have corresponding versions here, and function names should match the Haskell versions, so the linked Haskell documentation is a very useful reference to see what each function does. Even the function implementations themselves are almost identical to the Haskell versions (excluding changes based on C++ vs. Haskell language differences).

## Dependencies

This project relies on [GMP](https://gmplib.org/) and [Boost](https://www.boost.org/). It was tested with GMP 6.2.1 and Boost 1.80.

## Types

Aliases are defined in "types.h" for certain types to more closely match Haskell names (e.g. use `Maybe` instead of `std::optional`). Arbitrary precision integers and rationals are represented by GMP's `mpz_class` (aliased as `Integer`) and `mpq_class` (aliased as `Rational`), respectively. For high-precision real numbers, the `cpp_dec_float` class from Boost is used (aliased as `Real`), which is used to mimic the `FixedPrec` numbers in Haskell.

## Running the Program

Run `make main` to produce an executable called `main` in the "build" folder. Then, run `build/main <angle> <prec> [effort=25]`:

- `angle`: the program approximates a z-rotation gate for this angle (in radians).
- `prec`: the required precision in bits for the approximation.
- `effort`: the effort used for factoring. This is not required and defaults to 25.

These values correspond to the second, third, and fourth arguments of [gridsynth](https://hackage.haskell.org/package/newsynth-0.4.0.0/docs/Quantum-Synthesis-GridSynth.html#v:gridsynth) in Haskell.

For example, `build/main 1.2 10` will synthesize an approximate z-rotation gate of 1.2 radians to at least 10 bits of precision. The output will include the sequence of gates (as a string), along with the $\log_{0.5}$ of the approximation error. These values correspond to the gate sequence returned by [gridsynth_gates](https://hackage.haskell.org/package/newsynth-0.4.0.0/docs/Quantum-Synthesis-GridSynth.html#v:gridsynth_gates) and the error returned by [gridsynth_stats](https://hackage.haskell.org/package/newsynth-0.4.0.0/docs/Quantum-Synthesis-GridSynth.html#v:gridsynth_stats) in Haskell.

### Choosing Precision

The `prec` input controls how accurate the computed gate sequence has to be. Increasing precision increases runtime on a seemingly exponential basis for the C++ implementation. For example, on a 2020 Macbook Air, synthesis an angle of 1.7 radians took approximately:

- 4 seconds for 8 bits of precision
- 83 seconds for 9 bits
- 368 seconds for 10 bits
- 1288 seconds for 11 bits

This runtime will also vary for different angles, but using have too high of a precision value leads to very high runtimes, so starting with lower precisions and trying to work up is advisable. **This implementation is much slower than the Haskell version, which can still run in just a few seconds for 100s of bits of precision.** Some potential ways to improve the speed of this program are presented later.

### Setting Include Directory

The `INCLUDE_DIR` variable in the Makefile, which defaults to `/usr/local/include` may have to be changed depending on the system being used. This is where the GMP headers and the "boost" directory should be located. To set the variable, add `INCLUDE_DIR=<some path>` when invoking `make`, like `make main INCLUDE_DIR=/sample/path`.

### Setting Real Number Precision

The `Real` type defaults to having 100 digits of precision, to change this, set the `REAL_DIGITS` variable when building, like `make main REAL_DIGITS=50`. Note that building and running the tests like this may or may not cause errors due to small precision differences, since they were tested with 100 digits of precision.

**Note**: Based on experimentation, it seems like the number of digits of precision for `Real` (`cpp_dec_float`) refers to the total number of digits in a value (e.g. `123.456` has 6 total digits). On the other hand, the precision for the `FixedPrec` real numbers in Haskell seems to refer to the number of digits after the decimal point. This means that the two precision values aren't exactly the same. As a result, the Haskell heuristic for the number of digits to use, which appears in `gridsynth_stats`, may be an underestimate for C++.

## Running Tests

`make tests` will compile all of the tests, with outputs in the "build" folder that can be run individually. `make runtests` will compile and run all of the tests.

## Generating Documentation

Running `doxygen` (assuming it is installed) should generate HTML and Latex documentation (in "html" and "latex" directories, respectively). I haven't tried looking at the Latex documentation, but open up the "index.html" file in the "html" directory to view the HTML documentation.

## Areas of Future Improvement

This implementation is significantly slower than the Haskell implementation for larger precisions, and I had a few ideas on how to speed it up. In addition to these ideas, more testing should be done to identify which parts of the algorithm are taking up most of the runtime.

- Almost everything is passed by value, including lists, matrices, and lambda functions (as `std::function` values). If a long list or other object is passed around between functions, a lot of copying might be going on. Passing more things by const reference should should speed things up. This wasn't done initially because with all of the higher-order functions and lambdas in the code, I thought passing everything by value would make it easier to avoid bugs, but this is likely slowing execution down. Note that this is probably less of a problem for matrices, because almost all (if not all) of the matrices in the code are only 2x2.
  - A similar critique is true for functions like `utils::tail` that manipulate lists: if we take the tail of a very long list, almost all of the list is being copied, so some sort of range that points to the original list could be used instead to avoid copying.
- The numeric classes, like `RootTwo` or `Complex`, are also passed by value and stored as values almost all of the time (e.g. a list will be `List<QRootTwo>` instead of `List<const QRootTwo&>`, and the same for matrices). Since these numeric classes often contain arbitrary precision integers or rationals, which take up more space than normal primitive types, passing them around by reference might improve execution speed.
- Some recursive algorithms might be convertible to iterative ones, which could save space/time.
- Also refer to any `\todo` comments in the code.
