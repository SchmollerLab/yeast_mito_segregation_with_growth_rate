# +
#pip install aicsimageio
# -

import os
import numpy as np

import tifffile
from tqdm.notebook import tqdm

import pandas as pd
import dask

from IPython.core.debugger import set_trace

ND2_FILEPATH = '/volumes/roussou/Microscopy/CellAsic/Atp6_NG_x_Atp6_mKate_Deltapdr5_Chloramphenicol/Diploids_withCAP_washedafter6hrs_010623/raw_data/delta_pdr5_ng_x_mkate_capfor6hrs_010624.nd2'
EXP_FOLDERPATH = '/volumes/roussou/Microscopy/CellAsic/Atp6_NG_x_Atp6_mKate_Deltapdr5_Chloramphenicol/Diploids_withCAP_washedafter6hrs_010623/Positions/'
CHANNELS =  ['BF', 'BF1', 'mKate', 'NG']   # ['BF', 'BF1','None', 'NG', 'mKate','None'] #  ['BF', 'NG'] 

import nd2
nd2_data = nd2.ND2File(ND2_FILEPATH)
nd2_data.shape

nd2_dask = nd2_data.to_dask()

nd2_img = nd2_dask[0, 0, 0, 0]
nd2_img.shape

SizeT, SizeS, SizeZ, SizeC, SizeY, SizeX = nd2_data.shape
pos_pbar = tqdm(total=SizeS, desc='Position')
for s in range(SizeS):
#     if s > 0:
#         break
    images_path = os.path.join(EXP_FOLDERPATH, f'Position_{s+1}', 'Images')
    if not os.path.exists(images_path):
        continue
    tif_filenames = set()
    df_roi = None
    shifts = None
    for acdc_filename in os.listdir(images_path):
        acdc_filepath = os.path.join(images_path, acdc_filename)
        if acdc_filename.endswith('dataPrepROIs_coords.csv'):
            df_roi = pd.read_csv(acdc_filepath).set_index('description')
        elif acdc_filename.endswith('.tif'):
            tif_filenames.add(acdc_filename)
        elif acdc_filename.endswith('align_shift.npy'):
            shifts = np.load(acdc_filepath)
    basename = os.path.commonprefix(list(tif_filenames))
    if df_roi is not None:
        x0 = int(df_roi.at['x_left', 'value'])
        x1 = int(df_roi.at['x_right', 'value'])
        y0 = int(df_roi.at['y_top', 'value'])
        y1 = int(df_roi.at['y_bottom', 'value'])
    for c, channel in enumerate(CHANNELS):
        if channel == 'None':
            continue
        channel_data = nd2_dask[:, s, :, c]
        aligned_data = np.zeros(channel_data.shape, channel_data.dtype)
        frames_pbar = tqdm(total=len(channel_data), desc='Aligning', leave=False)
        for frame_i in range(len(shifts)):
            shifts_i = shifts[frame_i]
            img = channel_data[frame_i].compute()
            aligned_img = np.roll(img, tuple(shifts_i), axis=(1,2))
            aligned_data[frame_i] = aligned_img
            frames_pbar.update()
        frames_pbar.close()
        if df_roi is not None: 
            cropped_aligned_data = aligned_data[..., y0:y1, x0:x1] # .compute()
        else:
            cropped_aligned_data = aligned_data
        npz_filename = f'{basename}_{channel}_3D_aligned.npz'
        npz_filepath = os.path.join(images_path, npz_filename)
        tif_filename = f'{basename}_{channel}_3D.tif'
        tif_filepath = os.path.join(images_path, tif_filename)
        # np.savez_compressed(npz_filepath, cropped_aligned_data)
        tifffile.imwrite(tif_filepath, cropped_aligned_data)
        # set_trace()
    pos_pbar.update()

pos_pbar.close()
