#!/bin/bash
# Step 1: Generate gold-standard protein-protein proximity values
python calibrate_pairwise_distance.py --protein_file ./Examples/MuSIC_proteins.txt --outprefix ./Examples/output/test

# Step 2: Build random forest to predict protein-protein proximity from data embeddings
python random_forest_samples.py --outprefix ./Examples/output/test --protein_file ./Examples/MuSIC_proteins.txt --emd_files ./Examples/IF_image_embedding.csv ./Examples/APMS_embedding.MuSIC.csv --emd_label IF_emd APMS_emd --num_set 2 auto --n_samples 1000 

for ((fold = 1; fold <= 5; fold++))
do
    for ((IF_set = 1; IF_set <= 2; IF_set++))
    do
        python run_random_forest.py --outprefix ./Examples/output/test --fold $fold --emd_label IF_emd APMS_emd --train_set $IF_set 1 --n_jobs 60;
    done
done

python random_forest_output.py --outprefix ./Examples/output/test

# Step 3: Analyze proximity data to identify protein communities at progressive resolutions
cp ./Examples/MuSIC_predicted_proximity.txt ./Examples/output/test_predicted_proximity.txt
cp ./Examples/MuSIC_avgPred_ytrue.csv ./Examples/output/test_avgPred_ytrue.csv

python community_detection.py --outprefix ./Examples/output/test --path_to_clixo /cellar/users/y8qin/Modules/CliXO --clixo_i ./Examples/output/test_predicted_proximity.txt --clixo_a 0.01 --clixo_b 0.5 --clixo_m 0.008 --clixo_z 0.05 --min_diff 2 --path_to_alignOntology /cellar/users/y8qin/Modules/alignOntology-master --predict_nm_size --keep_all_files

bash ./Examples/output/test.sh