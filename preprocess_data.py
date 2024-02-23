import pandas as pd
import numpy as np
import h5py
import yaml
import argparse
import os
import progressbar


def load_config(config_path):
    try:
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config
    except Exception as e:
        print(f"Error loading configuration file: {e}")
        return {}


def main(config_path):
    config = load_config(config_path)
    print(f"Loaded Configuration:{config}")
    print('Splitting data by GO branch')
    
    if not os.path.exists(config['directories']['preprocessed_data']):
        os.makedirs(config['directories']['preprocessed_data'])
    
    path_to_train_set = config['directories']['raw_data_base'] + '/train/train_set.tsv'
    train_set_df = pd.read_csv(path_to_train_set, sep='\t')
    
    aspect_groups = train_set_df.groupby('aspect')
    # Create separate DataFrames for each 'aspect'
    aspect_dfs = {aspect: aspect_groups.get_group(aspect) for aspect in aspect_groups.groups}
    
    print("\nDataFrames split by 'aspect':")
    for aspect, aspect_df in aspect_dfs.items():
        if aspect == 'biological_process':
            BP_df = aspect_df
        elif aspect == 'cellular_component':
            CC_df = aspect_df
        else: #aspect == 'molecular_function':
            MF_df = aspect_df

    print(f"\nDataFrame for aspect Biological Processes with {len(BP_df)} elements")
    print(BP_df.head())

    print(f"\nDataFrame for aspect Cellular Component with {len(CC_df)} elements")
    print(CC_df.head())

    print(f"\nDataFrame for aspect Molecular Function with {len(MF_df)} elements")
    print(MF_df.head())
    
    path_to_train_embeddings = config['directories']['raw_data_base'] + '/train/train_embeddings.h5'
    embeddings_hf = h5py.File(path_to_train_embeddings, 'r')
    hf_dsets = list(embeddings_hf.keys())
    loa = []
    for key in hf_dsets:
        arr = np.array(embeddings_hf[key])
        loa.append(arr)
    train_embed_df = pd.DataFrame(loa)
    train_embed_df.insert(0, "Protein_ID", hf_dsets)
    train_embed_df.set_index("Protein_ID", inplace=True)
    
    print(train_embed_df.head())
    
    BP_Protein_IDs = BP_df.Protein_ID.unique()
    CC_Protein_IDs = CC_df.Protein_ID.unique()
    MF_Protein_IDs = MF_df.Protein_ID.unique()
    print('Number of unique proteins in BP: ', len(BP_Protein_IDs))
    print('Number of unique proteins in CC: ', len(CC_Protein_IDs))
    print('Number of unique proteins in MF: ', len(MF_Protein_IDs))
    
    train_embed_BP_df = train_embed_df.filter(items = BP_Protein_IDs, axis=0)
    train_embed_CC_df = train_embed_df.filter(items = CC_Protein_IDs, axis=0)
    train_embed_MF_df = train_embed_df.filter(items = MF_Protein_IDs, axis=0)
    
    train_embed_BP_df.to_pickle(config['directories']['preprocessed_data'] + '/train_embeddings_BiologicalProcesses.pkl')
    train_embed_CC_df.to_pickle(config['directories']['preprocessed_data'] + '/train_embeddings_CellularComponent.pkl')
    train_embed_MF_df.to_pickle(config['directories']['preprocessed_data'] + '/train_embeddings_MolecularFunction.pkl')
    
    train_protein_ids = np.loadtxt(config['directories']['raw_data_base'] + '/train/train_ids.txt', dtype = str)
    
    # Set the limit for label
    num_labels = 1500

    # Take value counts in descending order and fetch first 1500 `GO term ID` as labels
    labels = train_set_df['GO_term'].value_counts().index[:num_labels].tolist()

    # Fetch the train_terms data for the relevant labels only
    train_set_top1500terms = train_set_df.loc[train_set_df['GO_term'].isin(labels)]
    
    bar = progressbar.ProgressBar(maxval=num_labels, \
        widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    
    grouped = train_set_top1500terms.groupby('GO_term')['Protein_ID'].apply(set).to_dict()

    set_ids = [
        ('BiologicalProcesses', BP_Protein_IDs),
        ('CellularComponent', CC_Protein_IDs),
        ('MolecularFunction', MF_Protein_IDs),
    ]

    for title, ids in set_ids:

        id_size = ids.shape[0]
        label_set = np.zeros((id_size, num_labels))

        bar.start()
        for i, label in enumerate(labels):
            label_related_proteins = grouped.get(label, set())
            label_set[:, i] = [id in label_related_proteins for id in ids]
            bar.update(i + 1)
        bar.finish()

        labels_df = pd.DataFrame(data=label_set, columns=labels)
        labels_df.insert(0, "Protein_ID", ids)
        labels_df.set_index("Protein_ID", inplace=True)
        labels_df.to_pickle(f'{config["directories"]["preprocessed_data"]}/train_labels_{title}.pkl')

        print(title, labels_df.shape)
    
    datasets = ['BiologicalProcesses', 'MolecularFunction', 'CellularComponent']

    for dataset in datasets:
        label_df = pd.read_pickle(f'{config["directories"]["preprocessed_data"]}/train_labels_{dataset}.pkl')
        embedding_df = pd.read_pickle(f'{config["directories"]["preprocessed_data"]}/train_embeddings_{dataset}.pkl')

        sampled_embedding_df = embedding_df.sample(frac=0.75, random_state=5)
        sampled_label_df = label_df.loc[sampled_embedding_df.index]
        
        sampled_embedding_df.to_pickle(f'{config["directories"]["preprocessed_data"]}/75percent_train_embeddings_{dataset}.pkl')
        sampled_label_df.to_pickle(f'{config["directories"]["preprocessed_data"]}/75percent_train_labels_{dataset}.pkl')
        print(f'Saved 75% {dataset} dataset')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Data preprocessing")
    parser.add_argument('--config', type=str, default='./config.yaml', help='Path to the configuration file')
    
    args = parser.parse_args()

    main(args.config)

