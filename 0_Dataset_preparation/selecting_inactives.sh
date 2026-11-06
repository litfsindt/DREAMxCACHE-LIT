input="."
output="."
target_sample="50000"
seed="20260807"
workers="20"

python3 sample_random_opendel.py --input-dir ${input} --output-dir ${output} --target ${target_sample} --seed ${seed} --workers ${workers}

awk 'NR==FNR {if (NF) pgk[$1]=1; next} NF && !($1 in pgk)' PGK2_selection_7M.smi ${output}/opendel_random_${target_sample}.smi > inactives.smi
