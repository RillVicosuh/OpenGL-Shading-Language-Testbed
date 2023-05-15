# GLSL Testbed

This repository contains a simple testbed for your GLSL implementation. The purpose of this testbed is to help you test and verify your GLSL code by comparing your output with reference solutions.

## Usage

To run the testbed, use the following command:

```
./driver -i <input-file> [ -s <solution-file> ] [ -o <stats-file> ]
```

Where:

- `<input-file>` is a file with commands to run.
- `<solution-file>` (optional) is a file with a solution to compare with.
- `<stats-file>` (optional) is a file to dump statistics to rather than stdout.

For example, to run the test with input file `00.txt`:

```
./driver -i 00.txt
```

This will save the result to `output.png`. You may compare this result with a reference solution:

```
./driver -i 00.txt -s 00.png
```

This will compute and output a measure of the difference, such as "diff: 0.23"; to pass a test, this error must be below the test's threshold. An image is output to `diff.png`, which visually shows where the differences are in the results, helping you track down any discrepancies.

The `-o` flag is used for the grading script so that grading is not confused by debug print statements.
