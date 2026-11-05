dataset="dataset.csv"
task_type="classification"
output_dir="chemprop_results"
split_type="scaffold_balanced"
seed="0" # We use 0, 1, 2, 3, 4
message_hidden_dim="200"

chemprop train --data-path ${dataset} --task-type ${task_type} --output-dir ${output_dir} --split-type ${split_type} --data-seed ${seed} --message-hidden-dim ${message_hidden_dim}
