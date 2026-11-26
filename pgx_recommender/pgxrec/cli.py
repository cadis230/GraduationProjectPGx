"""Usage:
    python -m pgxrec.cli --cyp2c19 "*2/*17"
    python -m pgxrec.cli --cyp2c19 "*1/*1" --cyp2d6 "*4/*4"
"""
import argparse

from .recommend import format_result, recommend


def main() -> None:
    p = argparse.ArgumentParser(description="CPIC-based antidepressant ranking from CYP2C19 (+ optional CYP2D6) diplotypes")
    p.add_argument("--cyp2c19", required=True, help="diplotype, e.g. '*1/*17'")
    p.add_argument("--cyp2d6", help="optional diplotype, e.g. '*1/*4' or '*2/*2x2' (copy number notation as in CPIC table)")
    args = p.parse_args()
    print(format_result(recommend(args.cyp2c19, args.cyp2d6)))


if __name__ == "__main__":
    main()
