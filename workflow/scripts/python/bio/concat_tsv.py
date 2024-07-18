import pandas as pd
from common.tsv import write_tsv
from common.bed import sort_bed_numerically

# TODO this will use a boatload of memory unnecessarily. We can assume both
# inputs are sorted, and therefore we simply need to stream both of them and
# merge based on the first 3 columns (constant memory)

# use pandas here since it will more reliably account for headers
df = pd.concat([pd.read_table(i, header=0) for i in snakemake.input])  # type: ignore
write_tsv(snakemake.output[0], sort_bed_numerically(df, 3))  # type: ignore
