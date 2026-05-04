"""Command-line interface for both sequence alignment versions."""

import argparse
from pathlib import Path

from sequence_alignment.basic import align_basic
from sequence_alignment.efficient import align_efficient
from sequence_alignment.input_parser import read_problem
from sequence_alignment.metrics import run_with_metrics

ALIGNERS = {
    "basic": align_basic,
    "efficient": align_efficient,
}


def write_alignment_output(output_path, result, metrics):
    lines = [
        str(result.cost),
        result.aligned_x,
        result.aligned_y,
        "{:.6f}".format(metrics.elapsed_ms),
        "{:.6f}".format(metrics.peak_memory_kb),
    ]
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(mode, input_path, output_path):
    sequence_x, sequence_y = read_problem(input_path)
    aligner = ALIGNERS[mode]
    result, metrics = run_with_metrics(lambda: aligner(sequence_x, sequence_y))
    write_alignment_output(output_path, result, metrics)
    return result, metrics


def _build_parser(default_mode):
    parser = argparse.ArgumentParser(description="Run CSCI570 sequence alignment.")
    if default_mode is None:
        parser.add_argument("mode", choices=sorted(ALIGNERS))
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    return parser


def main(default_mode=None, argv=None):
    parser = _build_parser(default_mode)
    args = parser.parse_args(argv)
    mode = default_mode if default_mode is not None else args.mode
    run(mode, args.input_file, args.output_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
