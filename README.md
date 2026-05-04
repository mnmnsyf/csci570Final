# CSCI570 Final Project

This repository is for the CSCI570 sequence alignment final project. It includes the two required versions of the algorithm: the basic dynamic programming version and the memory-efficient version.

The current version is a runnable project skeleton. The parser and algorithm modules expose fixed interfaces so that each part can be implemented and checked independently before integration.

## Dependencies

Use Python standard library modules plus `psutil` for memory measurement. No other external Python libraries are required.

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

Before the parser and algorithms are completed, the program only verifies that the command-line flow can run and write an output file.

## Data Files

`data/SampleTestCases/` contains official sample inputs and expected outputs. Use these files for correctness checks while implementing the parser and algorithms.

`data/Datapoints/` contains official input files for runtime and memory experiments. Use these files to compare the basic and memory-efficient versions after both are implemented.

`docs/` is reserved for project documents such as summaries, reports, and assignment notes. Test inputs and experiment inputs should stay under `data/`.

## Testing Plan

Optional unit tests can be added under `tests/` as modules are implemented. The PDF-required correctness checks and experiments are listed in the Testing role below.

Run tests with:

```powershell
python -m unittest discover
```

## Group Split

核心实现分给三位同学，另设一位同学负责测试和实验验收。输入解析、基础算法和内存优化算法通过固定接口连接，完成后可以直接通过 `basic.py` 和 `efficient.py` 运行。

### 1. Input Parser

负责文件：`src/sequence_alignment/input_parser.py`

主要任务：

- 读取 assignment input file。
- 根据 base string 和 index 列表生成最终的两条输入字符串。
- 校验输入格式和 DNA 字符是否合法。
- 保持 penalty 不从 input file 读取；gap penalty 和 mismatch cost 已经固定在 `src/sequence_alignment/constants.py`。

需要实现的接口：

```python
def generate_sequence(base_sequence, insertion_indices):
  ...

def parse_problem_lines(lines):
  ...
```

最终 `read_problem(input_path)` 应该返回：

```python
(sequence_x, sequence_y)
```

验收方式：

- `read_problem(input_path)` 能读取 `data/SampleTestCases/` 和 `data/Datapoints/` 下的输入文件。
- 返回值必须是 `(sequence_x, sequence_y)`，两条字符串只能包含 `A`、`C`、`G`、`T`。
- 对 sample cases，可以把官方 output 的 aligned strings 去掉 `_` 后，和 parser 生成的两条字符串对比。
- Input Parser 不负责读取 penalty；gap penalty 和 mismatch cost 必须继续使用 `src/sequence_alignment/constants.py` 中的固定值。

### 2. Basic Algorithm

负责文件：`src/sequence_alignment/basic.py`

主要任务：

- 实现普通 dynamic programming sequence alignment。
- 构建完整 DP table。
- 根据 DP table backtrack 得到一组 optimal alignment。
- 使用 `GAP_PENALTY` 和 `MISMATCH_COSTS` 计算 cost。

需要实现的接口：

```python
def build_cost_table(sequence_x, sequence_y):
  ...

def align_basic(sequence_x, sequence_y):
  ...
```

`align_basic` 必须返回：

```python
AlignmentResult(cost, aligned_x, aligned_y)
```

验收方式：

- `align_basic(sequence_x, sequence_y)` 必须返回 `AlignmentResult(cost, aligned_x, aligned_y)`。
- `aligned_x` 和 `aligned_y` 长度必须相同。
- 两条 aligned strings 去掉 `_` 后，必须分别还原为 parser 输出的 `sequence_x` 和 `sequence_y`。
- 根据 aligned strings 重新计算一次 alignment cost，结果必须等于 `AlignmentResult.cost`。
- 跑 `data/SampleTestCases/input1.txt` 到 `input5.txt`，输出第一行 cost 必须和官方 `output1.txt` 到 `output5.txt` 的第一行一致。

### 3. Memory-Efficient Algorithm

负责文件：`src/sequence_alignment/efficient.py`

主要任务：

- 实现 memory-efficient sequence alignment。
- 使用 linear-space DP 计算分割点。
- 使用 divide-and-conquer/Hirschberg 思路恢复 alignment。
- 保证最终 cost 和 basic version 一致。

需要实现的接口：

```python
def linear_space_scores(sequence_x, sequence_y):
  ...

def align_efficient(sequence_x, sequence_y):
  ...
```

`align_efficient` 必须返回：

```python
AlignmentResult(cost, aligned_x, aligned_y)
```

验收方式：

- `align_efficient(sequence_x, sequence_y)` 必须返回和 basic version 相同结构的 `AlignmentResult`。
- 在所有 sample cases 上，efficient version 的 cost 必须和 basic version 的 cost 一致。
- 在 `data/Datapoints/` 上，efficient version 也应和 basic version 得到相同 cost。
- 对较大的 datapoints，记录 basic 和 efficient 的 runtime/memory，用于 summary/report 中比较内存优化效果。

### 4. Testing

负责目录：`tests/`

完成 correctness check、命令行输出检查，以及 runtime/memory 实验数据收集。

验收方式：

- 先确认 Input Parser、Basic Algorithm、Memory-Efficient Algorithm 各自的验收项已经通过。
- 检查命令行接口是否符合 PDF 要求：程序必须接收 2 个参数，分别是 input file path 和 output file path。
- 检查 output file 不存在但路径合法时，程序是否会自动创建 output file。
- 使用 `data/SampleTestCases/input1.txt` 到 `input5.txt` 运行 basic version，并和对应 `output1.txt` 到 `output5.txt` 对比。
- 对 sample outputs，重点检查输出第一行 minimum alignment cost 是否和官方 output 第一行一致。
- 检查输出文件是否符合 5 行格式：cost、aligned string 1、aligned string 2、time、memory。
- 检查 aligned string 1 和 aligned string 2 长度是否相同。
- 检查 aligned strings 去掉 `_` 后是否能还原 input parser 生成的原始字符串。
- 在 memory-efficient version 完成后，检查 efficient version 和 basic version 在同一个 input 上的 cost 是否一致。
- 使用 `data/Datapoints/` 跑 basic 和 efficient，收集 runtime 和 memory，用于 summary/report。
- 检查 gap penalty 和 mismatch penalties 使用的是 `src/sequence_alignment/constants.py` 中的固定值，而不是从 input file 读取。
- 三位实现同学提交代码后，统一运行完整检查，确认 parser、basic、efficient、CLI 可以连通。

可选辅助测试文件：

```text
tests/test_input_parser.py
tests/test_basic.py
tests/test_efficient.py
tests/test_cli.py
```

如果添加了 unit tests，可以用下面命令统一运行：

```powershell
python -m unittest discover
```

### Testing and Experiments

- 使用 `data/SampleTestCases/` 做 correctness check。
- 使用 `data/Datapoints/` 收集 basic 和 efficient 的 runtime/memory。
- 如需自动化实验，可以在 `scripts/` 下添加脚本，例如 `scripts/run_datapoints.py`。
- 三位实现同学完成接口后，Testing 同学统一通过下面命令检查整体流程：

```powershell
python basic.py data/SampleTestCases/input1.txt outputs/basic_input1.txt
python efficient.py data/SampleTestCases/input1.txt outputs/efficient_input1.txt
```

## Final Submission Checklist

- `basic.py` and `efficient.py` both run with exactly 2 command-line arguments.
- Sample test costs match the official outputs in `data/SampleTestCases/`.
- Output files use the required 5-line format.
- Runtime and memory results are collected on `data/Datapoints/`.
- Gap penalty and mismatch penalties remain hardcoded in `src/sequence_alignment/constants.py`.
