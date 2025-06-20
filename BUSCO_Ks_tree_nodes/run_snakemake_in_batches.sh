for i in {1..199}; do
    snakemake --config part=${i} out_of=200 --cores [number of cores] --scheduler greedy -q all --keep-going
done
