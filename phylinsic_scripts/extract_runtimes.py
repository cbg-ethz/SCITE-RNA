import argparse

def main():
    parser = argparse.ArgumentParser(description="Extract and combine runtimes in order.")
    parser.add_argument("--input_files", nargs="+", required=True, help="List of input files (already ordered).")
    parser.add_argument("--output_file", required=True, help="Output file to write runtimes.")
    args = parser.parse_args()

    # Write runtimes to output in order
    with open(args.output_file, "w") as out_f:
        for file in args.input_files:
            with open(file, "r") as in_f:
                runtime = in_f.read().strip()
                out_f.write(runtime + "\n")

if __name__ == "__main__":
    main()
