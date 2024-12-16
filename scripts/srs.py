from rpy2.robjects import pandas2ri
pandas2ri.activate()
import rpy2.robjects as robjects
import read_data as rd
import random as rand

import pandas as pd

srs_path = 'R\\SRS.R'

def srs_on_id_list(id_list: list, new_read_count: int, seed = 0) -> dict:
    """_summary_

    Args:
        id_list (list): _description_: List of read ids
        new_read_count (int): _description_: Count that SRS will normalize to
        seed (int, optional): _description_: Seed for SRS algorithm. Defaults to 0.

    Returns:
        dict: _description_: Dictionary of read ids with normalized counts
    """    
    read_count_dict = {}
    for id in id_list:
        if id not in read_count_dict:
            read_count_dict[id] = 0
        read_count_dict[id] += 1
    
    r = robjects.r
    r['source'](srs_path)
    srs_func_r = robjects.globalenv['SRS']
    
    # convert dict to R df
    df = pd.DataFrame.from_dict(read_count_dict, orient='index')
    df_r = pandas2ri.py2rpy(df)
    srs_result = srs_func_r(df_r, new_read_count, True, seed)
    df = pandas2ri.rpy2py(srs_result)

    srs_result_values = list(df.to_dict()['0'].values())
    read_count_dict_labels = list(read_count_dict.keys())
    normalized_read_count_dict = {}
    
    for idx, val in enumerate(read_count_dict_labels):
        normalized_read_count_dict[val] = srs_result_values[idx]
    
    normalized_id_list = []
    for id in normalized_read_count_dict.keys():
        current_count = normalized_read_count_dict[id]
        for i in range(int(current_count)):
            normalized_id_list.append(id)
            
    rand.shuffle(normalized_id_list)
        
    return normalized_id_list


def srs_on_read_data(data: rd.ReadData, new_read_count: int, seed = 0) -> rd.ReadData:
    """_summary_
    
    Perform SRS on the class hashtable of reads from a ReadData class instance.

    Args:
        data (rd.ReadData): _description_: Hash table of read data
        new_read_count (int): _description_: Count that SRS will normalize to
        seed (int, optional): _description_: Seed for SRS algorithm. Defaults to 0.

    Returns:
        rd.ReadData: _description_: New ReadData instance with normalized read counts
    """   
    
    if new_read_count >= len(data.reads):
        return data
    
    reads_processed = data.processed_read_count
    data = data.reads
    
    read_count_dict = {}
    hash_dict = {}
    for key in data.keys():
        current_read : rd.Read = data[key]
        current_assignment = current_read.get_assignment()
        if current_assignment not in read_count_dict:
            read_count_dict[current_assignment] = 0
        if current_assignment not in hash_dict:
            hash_dict[current_assignment] = []
        read_count_dict[current_assignment] += 1
        hash_dict[current_assignment].append(key)

    r = robjects.r
    r['source'](srs_path)
    srs_func_r = robjects.globalenv['SRS']
  
    df = pd.DataFrame.from_dict(read_count_dict, orient='index')
    df_r = pandas2ri.py2rpy(df)
    srs_result = srs_func_r(df_r, new_read_count, True, seed)
    df = pandas2ri.rpy2py(srs_result)

    srs_result_values = list(df.to_dict()['0'].values())
    read_count_dict_labels = list(read_count_dict.keys())
    normalized_read_count_dict = {}
    
    for idx, val in enumerate(read_count_dict_labels):
        normalized_read_count_dict[val] = srs_result_values[idx]
    
    normalized_read_data: rd.ReadData = rd.ReadData()
    for id in normalized_read_count_dict.keys():
        current_count = normalized_read_count_dict[id]
        current_hashes = hash_dict[id]
        li = len(current_hashes)
        for i in range(int(current_count)):
            random_hash = rand.choice(current_hashes)
            normalized_read_data.insert_hash_as_read(random_hash, id)
            current_hashes.remove(random_hash)
        hash_dict[id] = current_hashes
    
    normalized_read_data.processed_read_count = reads_processed
    return normalized_read_data