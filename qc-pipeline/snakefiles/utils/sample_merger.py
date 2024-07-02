# #drafting some methods for merging batches. Work in progress...

# def merge_bedsets(file_dict, batches, path):
#     duplicate_individuals = find_duplicate_individuals(file_dict["fam"], batches, path)
#     samples_to_remove = choose_samples_to_remove(duplicate_individuals, file_dict)
#     remove_samples(samples_to_remove)

# def find_duplicate_individuals(fam_files, batches, path):
#     dfs = []
#     for batch in batches:
#         fam_file = fam_files[batch]
#         df = pd.read_csv(fam_file, delim_whitespace=True, header=None, names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENOTYPE'])
#         df['Batch'] = batch
#         dfs.append(df)
#     combined_df = pd.concat(dfs, ignore_index=False)
#     duplicate_individuals = combined_df[combined_df.duplicated('IID', keep=False)]
#     duplicate_individuals.to_csv(path + '\duplicate_individuals.txt', sep='\t', index=False)
#     return duplicate_individuals


# def choose_samples_to_remove(duplicate_individuals, file_dict):
#     uniquie_duplicate_IIDs = duplicate_individuals.drop_duplicates(subset='IID')
#     remove_dict = {}
#     scores = []
#     for index, ind in uniquie_duplicate_IIDs.iterrows():
#         IID = ind['IID']
#         duplicates = duplicate_individuals[duplicate_individuals['IID'] == IID]
#         scores = []
#         batches = []
#         for index2, dup in duplicates:
#             batch = dup["batch"]
#             batches.append(batch)
#             fam_file = file_dict["fam"][batch]
#             score = keep_score(IID, fam_file)
#             scores.append(score)
#         max_index = scores.index(max(scores))
#         for i in len(batches):
#             if i != max_index:
#                 if batch in remove_dict:
#                     remove_dict[batch].append(IID)
#                 else:
#                     remove_dict[batch] = [IID]
#     return remove_dict

    
# def keep_score(IID, fam_file):
#     score = 0
#     family_in_batch = family_in_batch(IID, fam_file)
#     if family_in_batch:
#         score += 1
#     return score

# def family_in_batch(IID, fam_file):
#     fam_df = pd.read_csv(fam_file, delim_whitespace=True, header=None, names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENOTYPE'])
#     FID = fam_df[fam_df["IID"] == IID]["FID"].iat[0]
#     n_family_members = fam_df[fam_df["FID"]== FID].shape[0]
#     return n_family_members > 1

# # def remove_samples(samples_to_remove, path):
# #     for batch, IID in samples_to_remove.items():