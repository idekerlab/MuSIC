#!/bin/bash
echo "Starting toy example using 100 proteins with random embeddings."

# Step 1: Generate gold-standard protein-protein proximity values
python calibrate_pairwise_distance.py --protein_file ./Examples/toy/toy_proteins.txt --outprefix ./Examples/toy_output/toy

# Step 2: Build random forest to predict protein-protein proximity from data embeddings
python random_forest_samples.py --outprefix ./Examples/toy_output/toy --protein_file ./Examples/toy/toy_proteins.txt --emd_files ./Examples/toy/toy_IF_image_embedding.csv ./Examples/toy/toy_APMS_embedding.csv --emd_label IF_emd APMS_emd --n_samples 1000

for ((fold = 1; fold <= 5; fold++));
do
    python run_random_forest.py --outprefix ./Examples/toy_output/toy --fold $fold --emd_label IF_emd APMS_emd;
done

python random_forest_output.py --outprefix ./Examples/toy_output/toy

# Step 3: Analyze proximity data to identify protein communities at progressive resolutions
python community_detection.py --outprefix ./Examples/toy_output/toy --clixo_i ./Examples/toy_output/toy_predicted_proximity.txt --predict_nm_size --keep_all_files

bash ./Examples/toy_output/toy.sh