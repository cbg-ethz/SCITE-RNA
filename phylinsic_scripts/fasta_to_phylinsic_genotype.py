def main():
    infile = snakemake.input[0]
    outfile = snakemake.output[0]

    with open(infile, encoding="utf-8") as inf, open(outfile, "w", encoding="utf-8") as outf:
        first = True
        for line in inf:
            line = line.rstrip("\r\n")
            cols = line.split("\t")
            # If this is the header line (starts with "Mutation" or "mutation"),
            # skip it and continue.
            if first and cols[0].lower() == "mutation":
                first = False
                continue
            # Otherwise, drop the first column (mutation id) and write only genotypes
            # If the line has only one column, write an empty line
            if len(cols) <= 1:
                outf.write("\n")
            else:
                outf.write(" ".join(cols[1:]) + "\n")
            first = False

if __name__ == "__main__":
    main()