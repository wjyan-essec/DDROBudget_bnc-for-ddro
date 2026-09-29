"""Run BB-CONN or BC-INT-CONN on a converted TSP MPS/AUX instance."""

import argparse
import csv
from pathlib import Path

from globals import Config, Tracker
from run_bnc import read_instance_metadata
from tsp_problem_class import TSPBnCProblem, TSP_RESULT_FIELDS


REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_CSV = REPOSITORY_ROOT / "results" / "tsp_summary.csv"


def parse_command_line_arguments():
    parser = argparse.ArgumentParser(
        description="Run TSP DDRO BnC with mandatory connectivity cuts."
    )
    parser.add_argument("--instance_file", required=True)
    parser.add_argument(
        "--cuts", required=True, choices=["branchandbound", "intersection"]
    )
    parser.add_argument("--projected", type=int, default=1, choices=[0, 1])
    parser.add_argument(
        "--separation",
        default="integer",
        choices=["integer", "fractional"],
    )
    parser.add_argument("--time_lim", type=float, default=60.0)
    parser.add_argument("--output_csv", default=str(DEFAULT_OUTPUT_CSV))
    parser.add_argument("--verbose_level", type=int, default=0, choices=[0, 1, 2, 3])
    parser.add_argument("--max_cuts", type=int, default=20)
    parser.add_argument("--cplex_cuts", type=int, default=-1, choices=[-1, 0])
    parser.add_argument("--only_root_node", type=int, default=0, choices=[0, 1])
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument("--write_lps", type=int, default=0, choices=[0, 1])
    return parser.parse_args()


def validate_output_csv(output_csv):
    output_path = Path(output_csv).resolve()
    if not output_path.exists() or output_path.stat().st_size == 0:
        return
    with output_path.open(newline="", encoding="utf-8") as stream:
        header = next(csv.reader(stream), [])
    if header != list(TSP_RESULT_FIELDS):
        raise ValueError(
            f"Existing TSP result CSV has an incompatible header: {output_path}."
        )


def write_result_csv(output_csv, result):
    output_path = Path(output_csv).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not output_path.exists() or output_path.stat().st_size == 0
    with output_path.open("a", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=TSP_RESULT_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow(result)
    return output_path


def main():
    args = parse_command_line_arguments()
    metadata = read_instance_metadata(args.instance_file)
    metadata["problem"] = "TSP"
    source_file = metadata.get("source_file")
    if not source_file:
        raise ValueError(
            "The TSP MPS must have a conversion_manifest.csv entry with source_file. "
            "Generate it with convert_sp_to_bobilib.py."
        )
    source_path = Path(source_file)
    if not source_path.is_absolute():
        source_path = REPOSITORY_ROOT / source_path
    if not source_path.is_file():
        raise FileNotFoundError(f"Original SP source file not found: {source_path}")
    bound_time = metadata["bound_time"]
    model_time_limit = max(0.0, args.time_lim - bound_time)

    intersection = args.cuts == "intersection"
    config = Config(
        instance_type="bobilib",
        instance_file=args.instance_file,
        lower_level="general",
        projected=bool(args.projected),
        separation=args.separation,
        time_lim=model_time_limit,
        bound_time=bound_time,
        total_time_lim=args.time_lim,
        instance_metadata=metadata,
        verbose_level=args.verbose_level,
        max_cuts=args.max_cuts,
        only_root_node=bool(args.only_root_node),
        cplex_cuts=args.cplex_cuts,
        tolerance=args.tolerance,
        write_lps=bool(args.write_lps),
        interdiction_cuts=False,
        intersection_cuts=intersection,
        nogood_cuts=False,
        branchandbound=not intersection,
        tsp_source_file=str(source_path),
    )
    tracker = Tracker()
    validate_output_csv(args.output_csv)

    print("\n%%%%%%%%%% TSP parameter information %%%%%%%%%%")
    print(f"Instance file: {args.instance_file}")
    print(f"Original SP file: {source_path}")
    print(f"Method: {'BC-INT-CONN' if intersection else 'BB-CONN'}")
    print(f"Total time limit: {args.time_lim} seconds")
    if intersection:
        print(f"Intersection max cuts per node: {args.max_cuts}")
    print("Connectivity cut limit: none")
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

    problem = TSPBnCProblem(config, tracker)
    result = problem.solve()
    output_path = write_result_csv(args.output_csv, result)
    print(f"Result CSV: {output_path}")


if __name__ == "__main__":
    main()
