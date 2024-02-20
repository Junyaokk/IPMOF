import pandas as pd
import numpy as np
import datetime
import warnings
import inspect
warnings.filterwarnings("ignore")


def load_data(data_dir):
    file_names = ['sku_info.csv', 'region_info.csv', 'ord_info.csv', 'trans_info.csv']
    dfs = {file_name.split('.')[0]: pd.read_csv(data_dir + file_name) for file_name in file_names}

    ord_info_df = dfs['ord_info']
    ord_info_df['ord_date'] = pd.to_datetime(ord_info_df['ord_date']).dt.date
    ord_info_df = ord_info_df.sort_values(by='ord_date').reset_index(drop=True)

    return dfs['sku_info'], dfs['region_info'], ord_info_df, dfs['trans_info']


def get_sku_info(sku_info_df):
    sku_set = sorted(sku_info_df['sku_id'].tolist())
    sku_weight_dict_dict = sku_info_df.set_index('sku_id')['sku_weight'].to_dict()
    return sku_set, sku_weight_dict_dict


def get_supply_info(region_info_df):
    rdc_id = region_info_df.at[0, 'rdc_id']
    fdc_id = region_info_df.at[0, 'fdc_id']
    region_set = region_info_df['region_id'].unique().tolist()
    return rdc_id, fdc_id, region_set


def get_trans_info(trans_info_df, rdc_id, fdc_id):
    first_trans_df = trans_info_df[(trans_info_df['type'] == 'transshipment') & (trans_info_df['to'] == fdc_id)]
    second_trans_df = trans_info_df[(trans_info_df['type'] == 'fulfillment') & (trans_info_df['from'].isin([rdc_id, fdc_id]))]

    first_trans_dict = {rdc_id: {fdc_id: {'coef': first_trans_df['coef'].values[0], 
                                          'intercept': first_trans_df['intercept'].values[0]}}}
    second_trans_dict = {dc: {row['to']: {'coef': row['coef'], 
                                          'intercept': row['intercept']} 
                                          for _, row in group.iterrows()} 
                         for dc, group in second_trans_df.groupby('from')}

    return first_trans_dict, second_trans_dict


def count_unique_orders(data, key):
    unique_set = set()
    index_dict = data.get(key, {})
    
    for product_dict in index_dict.values():
        for order_list in product_dict.values():
            unique_set.update(order_list)
    
    return len(unique_set)


def remove_empty_dict_keys(dct):
    for k, v in list(dct.items()):
        if isinstance(v, dict):
            remove_empty_dict_keys(v)
            if not bool(v):
                del dct[k]
        elif not v:
            del dct[k]

def fill_missing_dates(sales_dict):
    min_date = min(sales_dict.keys())
    max_date = max(sales_dict.keys())
    
    filled_dict = {}
    
    current_date = min_date
    while current_date <= max_date:
        if current_date in sales_dict:
            filled_dict[current_date] = sales_dict[current_date]
        else:
            filled_dict[current_date] = 0
        
        current_date += datetime.timedelta(days=1)
    
    return filled_dict


def generate_new_data_by_sample(data, sample_size, num_samples):
    new_data = []
    
    for _ in range(num_samples):
        samples = np.random.choice(data, size=sample_size, replace=True)
        num = np.sum(samples)
        new_data.append(num)
    
    return new_data


def generate_new_data_by_sequence(data, num_elements):
    new_data = []
    
    for i in range(len(data) - num_elements + 1):
        subset = data[i:i + num_elements]
        sorted_subset = sorted(subset)
        num = sum(sorted_subset)
        new_data.append(num)
    
    return new_data


def closest_negative_to_zero(lst):
    closest_negative = None
    for num in lst:
        if num < 0:
            if closest_negative is None or abs(num) < abs(closest_negative):
                closest_negative = num
    return closest_negative


