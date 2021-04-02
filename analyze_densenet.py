import os
import pickle
from random import randint
from df_utils import *
from file_utils import *
import random
import argparse


parser = argparse.ArgumentParser(description='Split images into folds and analyze densenet embedding file.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--img_info', help='Full path to image embedding info file.')
parser.add_argument('--img_emd', help='Full path to image embedding file.')
parser.add_argument('--prefix', help='densenet embedding prefix E.g. densenet_raw_1024D')
args = parser.parse_args()

workdir = args.outprefix
img_info = args.img_info
img_emd = args.img_emd
label = args.prefix

# Split into six folds
print(img_info.shape)
if not os.path.exists('{}/densenet_emb_sim'.format(workdir)):
    os.mkdir('{}/densenet_emb_sim'.format(workdir))
for i in range(1, 7):
    if not os.path.exists('{}/densenet_emb_sim/fold_{}'.format(workdir, i)):
        os.mkdir('{}/densenet_emb_sim/fold_{}'.format(workdir, i))

genelist = list(set(img_info['gene_names']))
seen_img = dict()
for fold in range(1, 7):
    imgid_list = []
    for g in genelist:
        if fold is 1:
            tmp_imgid = img_info[img_info['gene_names'] == g].index.values
            selected_id = tmp_imgid[randint(0, len(tmp_imgid) - 1)]
            imgid_list.append(selected_id)
            seen_img[g] = list()
            seen_img[g].append(selected_id)
        else:
            tmp_imgid = img_info[img_info['gene_names'] == g].index.values
            tmp_seen = seen_img[g]
            tmp_unseen = list(set(tmp_imgid) - set(tmp_seen))
            # If all images are seen, sample from all images
            if len(tmp_unseen) == 0:
                imgid_list.append(tmp_seen[randint(0, len(tmp_seen) - 1)])
            # If there are unseen images, sample from unseen ones
            else:
                selected_id = tmp_unseen[randint(0, len(tmp_unseen) - 1)]
                imgid_list.append(selected_id)
                seen_img[g].append(selected_id)
    perFold_gname_to_imgID = dict(zip(genelist, imgid_list))
    perFold_imgID_to_gname = dict(zip(imgid_list, genelist))
    save_obj(perFold_gname_to_imgID, '{}/densenet_emb_sim/fold_{}/gname_to_imgID.dict.pkl'.format(workdir, fold))
    save_obj(perFold_imgID_to_gname, '{}/densenet_emb_sim/fold_{}/imgID_to_gname.dict.pkl'.format(workdir, fold))

gname_to_imgID = load_obj('{}/densenet_emb_sim/fold_1/gname_to_imgID.dict.pkl'.format(workdir))
sample_img_idx = []
for g in genelist:
    tmp_imgid = img_info[img_info['gene_names'] == g].index.values
    gsample = []
    if len(tmp_imgid) == 2:
        for tidx in tmp_imgid:
            if tidx == gname_to_imgID[g]:
                gsample += [tidx] * 2
            else:
                gsample += [tidx] * 3
    elif len(tmp_imgid) == 3:
        for tidx in tmp_imgid:
            if tidx == gname_to_imgID[g]:
                gsample += [tidx] * 1
            else:
                gsample += [tidx] * 2
    elif len(tmp_imgid) == 4:
        rand_id = np.random.choice(tmp_imgid)
        gsample = list(tmp_imgid)
        gsample.append(rand_id)
    else:
        for tidx in tmp_imgid:
            if tidx != gname_to_imgID[g]:
                gsample.append(tidx)
    random.shuffle(gsample)
    sample_img_idx.append(gsample)

for idx in range(5):
    fold = idx + 2
    fold_imgid = [x[idx] for x in sample_img_idx]
    perFold_gname_to_imgID = dict(zip(genelist, fold_imgid))
    perFold_imgID_to_gname = dict(zip(fold_imgid, genelist))
    save_obj(perFold_gname_to_imgID, '{}/densenet_emb_sim/fold_{}/gname_to_imgID.dict.pkl'.format(workdir, fold))
    save_obj(perFold_imgID_to_gname, '{}/densenet_emb_sim/fold_{}/imgID_to_gname.dict.pkl'.format(workdir, fold))

# determine similarities for each fold (cosine, Manhattan, Pearson, Spearman, Kendall, Euclidean)
for fold in range(1, 7):
    output_dir = '{}/densenet_emb_sim/fold_{}'.format(workdir, fold)
    imgID_to_gname = load_obj('{}/densenet_emb_sim/fold_{}/imgID_to_gname.dict.pkl'.format(workdir, fold))
    imgID = [x for x in imgID_to_gname]
    gname = [imgID_to_gname[x] for x in imgID_to_gname]
    df = load_obj(img_emd).loc[imgID]
    df.index = gname

    cos_sim = cosine_similarity_scaled(df)
    save_obj(cos_sim, '{}/{}.cosine.scaled.pkl'.format(output_dir, label))
    mht_sim = manhattan_similarity(df)
    save_obj(mht_sim, '{}/{}.manhattan.scaled.pkl'.format(output_dir, label))
    euc_sim = euclidean_similarity(df)
    save_obj(euc_sim, '{}/{}.euclidean.scaled.pkl'.format(output_dir, label))
    pcorr = pearson_scaled(df)
    save_obj(pcorr, '{}/{}.pearson.scaled.pkl'.format(output_dir, label))
    scorr = spearman_scaled(df)
    save_obj(scorr, '{}/{}.spearman.scaled.pkl'.format(output_dir, label))
    kcorr = kendall_scaled(df)
    save_obj(kcorr, '{}/{}.kendall.scaled.pkl'.format(output_dir, label))