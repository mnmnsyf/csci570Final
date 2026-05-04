# CSCI570 Final Project

This repository contains the Python implementation framework for the CSCI570 sequence alignment final project. The project includes two versions of the algorithm discussed in class: the basic dynamic programming version and the memory-efficient version.

The current codebase is organized as a runnable scaffold for group development. The parser and algorithm modules contain TODO entry points so that each part can be implemented and tested independently.

## Dependencies

The implementation should use Python standard library modules plus `psutil` for memory measurement. No other external Python libraries are required.

Install dependencies with:

```powershell
pip install -r requirements.txt
```

## Project Layout

```text
basic.py                         # command-line entry point for the basic version
efficient.py                     # command-line entry point for the memory-efficient version
src/sequence_alignment/
  alignment.py                   # shared AlignmentResult data model
  basic.py                       # basic dynamic programming algorithm
  cli.py                         # shared command-line interface and output writer
  constants.py                   # gap penalty and mismatch cost table
  efficient.py                   # memory-efficient alignment algorithm
  input_parser.py                # input parser and sequence generator
  metrics.py                     # runtime and memory measurement helpers
scripts/
  README.md                      # helper scripts for experiments and result checks
tests/
  README.md                      # unit tests and regression tests
data/
  SampleTestCases/               # official inputs with expected outputs
  Datapoints/                    # official inputs for performance experiments
docs/
  Summary.docx                   # project summary or report materials
```

## Running the Program

Both top-level entry points use the assignment-style command format:

```powershell
python basic.py <input_file> <output_file>
python efficient.py <input_file> <output_file>
```

Example:

```powershell
python basic.py data/SampleTestCases/input1.txt outputs/basic_input1.txt
python efficient.py data/SampleTestCases/input1.txt outputs/efficient_input1.txt
```

The output file should contain five lines:

```text
alignment cost
aligned first sequence
aligned second sequence
time in milliseconds
memory in KB
```

Until the TODO modules are completed, the program only verifies that the command-line flow can run and write an output file.

## Data Files

`data/SampleTestCases/` contains official sample inputs and expected outputs. These files should be used for correctness checks while implementing the parser and algorithms.

`data/Datapoints/` contains official input files for runtime and memory experiments. These files are also useful for comparing the basic and memory-efficient versions after both are implemented.

`docs/` is reserved for project documents such as summaries, reports, and assignment notes. Test inputs and experiment inputs should stay under `data/`.

## Testing Plan

Add tests under `tests/` as modules are implemented. Recommended checks include:

- input parsing and sequence generation
- alignment cost calculation
- output file format
- basic version cost compared with sample outputs
- efficient version cost compared with the basic version
- runtime and memory collection on files under `data/Datapoints/`

Run tests with:

```powershell
python -m unittest discover
```

## Suggested Group Split

- Input parsing: implement `src/sequence_alignment/input_parser.py`.
- Basic algorithm: implement `src/sequence_alignment/basic.py`.
- Memory-efficient algorithm: implement `src/sequence_alignment/efficient.py`.
- Testing: add unit tests under `tests/` and compare results with `data/SampleTestCases/`.
- Experiments and report: add scripts under `scripts/` for collecting results from `data/Datapoints/`.

If the course autograder requires specific filenames such as `basic_*.py` and `efficient_*.py`, keep the package structure unchanged and copy or rename the two top-level entry files as needed.
