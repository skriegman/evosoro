from glob import glob
import pandas as pd
from evosoro.tools.data_analysis import get_all_data, combine_experiments, plot_time_series

# directory names
exp_name_1 = "evo"
exp_name_2 = "evo_devo"

exp1_files = glob('./{}/run*/bestSoFar/bestOfGen.txt'.format(exp_name_1))
exp2_files = glob('./{}/run*/bestSoFar/bestOfGen.txt'.format(exp_name_2))

combined_df = combine_experiments([get_all_data(exp1_files), get_all_data(exp2_files)],
                                  [exp_name_1, exp_name_2])

combined_df.to_csv('{0}_vs_{1}_results.csv'.format(exp_name_1, exp_name_2))

combined_df = pd.read_csv('{0}_vs_{1}_results.csv'.format(exp_name_1, exp_name_2))

plot_time_series(combined_df, "evo_vs_evo_devo")
