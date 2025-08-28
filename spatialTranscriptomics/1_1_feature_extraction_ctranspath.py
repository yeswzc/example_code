##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## @ Danh-Tai HOANG
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import pickle
from loguru import logger
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os,sys,time,platform

UTIL_PATH = '/data/wuz6/project/25.pedHGG_spatial/06.Deeplearning/Tai_DeepPT_11125591/11slide_processing/'
sys.path.append(UTIL_PATH)
from utils_preprocessing import init_random_seed
from utils_preprocessing import *
#import openslide
from PIL import Image

import torch, torchvision
from torchvision import transforms
from torchvision.models import resnet50
import torch.nn as nn
from torch.utils.data import Dataset
import timm
#from ctran import ctranspath #uses modified timm need to set a different ENV

##        
## check available device
device = (torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu'))
logger.info(f"device:, {device}")
#os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

init_random_seed(random_seed=42)
resnet50_file = "/data/wuz6/project/25.pedHGG_spatial/06.Deeplearning/Tai_DeepPT_11125591/ResNet50_IMAGENET1K_V2.pt"
ctranspath_file = "/data/wuz6/project/25.pedHGG_spatial/06.Deeplearning/foundataionModels/ctranspath.pth"
##======================================================================================================
#path2storage = "/data/wuz6/project/25.pedHGG_spatial/06.Deeplearning/" #sys.argv[1]
#project = "pHGG" #sys.argv[2]
i_slide = int(sys.argv[1])
#print(f"project: {project}, i_slide: {i_slide}")
logger.info(f"Sample {i_slide}")


extract_raw_tile = True
save_tile_file = True
plot_mask = True

extract_CTransPath_features = True  #different conda (timm v0.5.4 modified)

extract_resnet50_features = False 
extract_MUSK_features = False #differnt conda
extract_UNI_features = False #large memory


if extract_CTransPath_features:
    ctranspath_PATH = '/data/wuz6/project/25.pedHGG_spatial/06.Deeplearning/foundataionModels/TransPath/'
    sys.path.append(ctranspath_PATH)
    from ctran import ctranspath


if extract_MUSK_features:
    #if(device != "cuda"):
    #    logger.info("ERROR: MUSK requires cuda GPU. Exit....")
    #    exit()
    from musk import utils, modeling
    from timm.models import create_model
    import torchvision
    from PIL import Image
    from timm.data.constants import IMAGENET_INCEPTION_MEAN, IMAGENET_INCEPTION_STD
    from huggingface_hub import login
    with open('hugging_face_token', 'r') as file:
    # Read the first line and strip whitespace
        mytoken = file.readline().strip()
    logger.info(f"my token = why would I tell you?")
    login(token = mytoken)

evaluate_edge = True
#evaluate_color = False

mag_assumed = 20
mag_selected = 20

tile_size = 2*int(225/2)
tile_size_half = int(tile_size/2)

mask_downsampling = 16
mask_tile_size = int(np.ceil(tile_size/mask_downsampling))
logger.info(f"tile size: {tile_size }mask_tile_size: {mask_tile_size}")

##---------------------------------------
#path2slide = path2storage #+ f"{project}_slides_data/"
#print("path2slide:", path2slide)
#path2meta = "../10metadata/"

path2output = "./"    
#path2ST = "/data/wuz6/project/25.pedHGG_spatial/06.Deeplearning/01.rna_preprocess/"
path2ST = f"/data/wuz6/data/aldape_lab/spatial10xg/spaceranger_v201_he/"
path2mask = f"{path2output}/mask/"
path2features = f"{path2output}/features/"
path2selected_barcode = f"{path2output}/selected_barcode/"

#metadata = pd.read_csv(f"../10match_slide_rna/{project}_slide_matched.csv")
#metadata = pd.read_csv(f"{path2meta}{project}_slide_matched.csv")
#metadata_file = "/data/wuz6/project/25.pedHGG_spatial/06.Deeplearning/01.rna_preprocess/00.Data_locations_v20250122.txt"
#metadata_file = f"/data/wuz6/data/aldape_lab/spatial10xg/spaceranger_v201_he/Data_locations_v20250516.txt"
metadata_file = f"/data/wuz6/data/aldape_lab/spatial10xg/spaceranger_v201_he/Data_locations_v20250716.txt.withImage"
metadata = pd.read_table(metadata_file)
#get the last three rows

slide_names = metadata.Outputname.values
slide_file_names = metadata.HE.values
start_time = time.time()

slide_name = slide_names[i_slide]
slide_file_name = slide_file_names[i_slide]

#features/AD01_ctranspath.npy
#skip if feature file already exists
if os.path.exists(f"{path2features}{slide_name}_ctranspath.npy"):
    logger.info(f"Feature file {path2features}{slide_name}_ctranspath.npy already exists. Skipping sample {slide_name}.")
    logger.info(f"finished -- i_slide: {i_slide}, total time: {int(time.time() - start_time)}")
    exit()
#skip if image file does not exist
if not os.path.exists(str(slide_file_name)):
    logger.info(f"Image file {slide_file_name} does not exist. Skipping sample {slide_name}.")
    logger.info(f"finished -- i_slide: {i_slide}, total time: {int(time.time() - start_time)}")
    sys.exit(1)  # Exit with error

logger.info(f"slide_name: {slide_name}. slide_file_name: {slide_file_name},")

#head /data/wuz6/data/aldape_lab/spatial10xg/spaceranger_v201_he/AD01/outs/spatial/tissue_positions.csv
#barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres
#GTCACTTCCTTCTAGA-1,0,0,0,1584,2050
#TAGATGCGTGGTTGAA-1,0,1,1,1776,2178

tissue_potision_file = f"{path2ST}/{slide_name}/outs/spatial/tissue_positions.csv"
tissue_potision = pd.read_csv(tissue_potision_file)
#only keep in_tissue == 1
tissue_potision = tissue_potision[tissue_potision.in_tissue == 1]
tissue_potision = tissue_potision.reset_index(drop=True)
logger.info(f"tissue_potision.shape: {tissue_potision.shape}")




##--------------------------
if not os.path.exists(path2mask):
    os.makedirs(path2mask)

if not os.path.exists(path2features):
    os.makedirs(path2features)

if not os.path.exists(path2selected_barcode):
    os.makedirs(path2selected_barcode)
##======================================================================================================


if save_tile_file:
    ## create tile_folder:
    tile_folder = f"{path2output}/tiles/{slide_name}/"
    logger.info(f"tile_folder: {tile_folder}")
    if not os.path.exists(tile_folder):
        os.makedirs(tile_folder)

##====================================================================================================== 
#path2slide = "/data/wuz6/project/25.pedHGG_spatial/06.Deeplearning/"
#print(f"{path2slide}{slide_file_name}")
#slide = openslide.OpenSlide(f"{path2slide}{slide_file_name}")


Image.MAX_IMAGE_PIXELS = None  # Disable the limit completely
im = Image.open(slide_file_name)

## magnification max
#if openslide.PROPERTY_NAME_OBJECTIVE_POWER in slide.properties:
#    mag_max = slide.properties[openslide.PROPERTY_NAME_OBJECTIVE_POWER]
#    print("mag_max:", mag_max)
#    mag_original = mag_max
#else:
#    print("[WARNING] mag not found, assuming: {mag_assumed}")
#    mag_max = mag_assumed
#    mag_original = 0
mag_max = mag_assumed
mag_original = mag_max

## downsample_level
downsampling = int(int(mag_max)/mag_selected)
logger.info(f"downsampling: {downsampling}")

##====================================================================================================== 
## slide size at largest level (level=0)
#px0, py0 = slide.level_dimensions[0]
px0, py0 = im.size
tile_size0 = int(tile_size*downsampling)
logger.info(f"px0: {px0}, py0: {py0}, tile_size0: {tile_size0}")

n_rows,n_cols = 1+int(np.ceil(py0/tile_size0)), 1+int(np.ceil(px0/tile_size0))
logger.info(f"n_rows: {n_rows}, n_cols: {n_cols}")
logger.info(f"spot tissue max pxl_row_in_fullres: {max(tissue_potision.pxl_row_in_fullres)}, max pxl_col_in_fullres: {max(tissue_potision.pxl_col_in_fullres)}")
logger.info(f"spot tissue min pxl_row_in_fullres: {min(tissue_potision.pxl_row_in_fullres)}, min pxl_col_in_fullres: {min(tissue_potision.pxl_col_in_fullres)}")


#logger.info(f"n_tiles_total: {n_tiles_total}")
#set min_pxl_row to min(tissue_potision.pxl_row_in_fullres) or 0
min_pxl_row = min(tissue_potision.pxl_row_in_fullres) if min(tissue_potision.pxl_row_in_fullres) < 0 else 0
min_pxl_col = min(tissue_potision.pxl_col_in_fullres) if min(tissue_potision.pxl_col_in_fullres) < 0 else 0

max_spot_col = int((max(tissue_potision.pxl_col_in_fullres)-min_pxl_col)/tile_size0)
max_spot_row = int((max(tissue_potision.pxl_row_in_fullres)-min_pxl_row)/tile_size0)

if(max_spot_col > n_cols):
     logger.info(f"Error: max spot col = {max_spot_col}, but max image row = {n_cols}", file=sys.stderr)
     n_cols = max_spot_col
     exit(0)
if(max_spot_row > n_rows):
     logger.info(f"Warning: max spot row = {max_spot_row}, but max image row = {n_rows}", file=sys.stderr)
     n_rows = max_spot_row     
#     exit(0)

n_tiles_total = n_rows*n_cols 
##====================================================================================================== 
img_mask = np.full((int(np.ceil((n_rows)*mask_tile_size)),int(np.ceil((n_cols)*mask_tile_size)),3),255).astype(np.uint8)
mask =     np.full((int(np.ceil((n_rows)*mask_tile_size)),int(np.ceil((n_cols)*mask_tile_size)),3),255).astype(np.uint8)

i_tile = 0
tiles_list = []
tiles_files_list = []
k_barcode_selected = []

if extract_raw_tile:
    import utils_color_norm
    color_norm = utils_color_norm.macenko_normalizer()
    ## evaluate tile using color_norm
    edge_mag_thrsh = 15  
    edge_fraction_thrsh = 0.5
    for k in range(tissue_potision.shape[0]):
        #logger.info(f"row: {row}/{n_rows}")
        barcode = tissue_potision.barcode[k]
        array_row = tissue_potision.array_col[k]
        array_col = tissue_potision.array_row[k]
        pxl_row_in_fullres = tissue_potision.pxl_row_in_fullres[k]
        pxl_col_in_fullres = tissue_potision.pxl_col_in_fullres[k]

        box = (pxl_col_in_fullres-tile_size_half, (pxl_row_in_fullres-tile_size_half), pxl_col_in_fullres+tile_size_half, (pxl_row_in_fullres+tile_size_half))
        tile = im.crop(box)
        ###mask is down sampling image to smaller magnification
        #row, col = int(pxl_row_in_fullres/tile_size0), int(pxl_col_in_fullres/tile_size0)
        row = int((pxl_row_in_fullres - min_pxl_row) / tile_size0)
        col = int((pxl_col_in_fullres - min_pxl_col) / tile_size0)
        logger.info(f"k:{k}:row:{row}:col:{col}:{barcode}:array_row/col:{array_row}:{array_col}:{pxl_row_in_fullres}:{pxl_col_in_fullres}")
        
                # Skip if row or col is negative or out of bounds
        #if row < 0 or col < 0 or row >= n_rows or col >= n_cols:
        #    logger.info(f"Skipping tile at k={k} due to invalid row/col: row={row}, col={col}")
        #    continue
        #if mask_tile_size == 0:
        #    logger.info(f"Skipping tile at k={k} due to mask_tile_size=0")
        #    continue

        mask_tile = np.array(tile.resize((mask_tile_size, mask_tile_size)))
        img_mask[int(row*mask_tile_size):int((row+1)*mask_tile_size),\
                 int(col*mask_tile_size):int((col+1)*mask_tile_size),:] = mask_tile
    
        ###
        tile = np.array(tile)
        ## evaluate tile
        select = evaluate_tile(tile, edge_mag_thrsh, edge_fraction_thrsh)
        if select == 1:
            k_barcode_selected.append(k)
            ## 2022.09.08: color normalization:
            tile_norm = Image.fromarray(color_norm.transform(tile))
            #tile_norm = Image.fromarray(tile) #skip norm
        
            mask_tile_norm = np.array(tile_norm.resize((mask_tile_size, mask_tile_size)))
      
            mask[int(row*mask_tile_size):int((row+1)*mask_tile_size),\
                 int(col*mask_tile_size):int((col+1)*mask_tile_size),:] = mask_tile_norm
        
            tiles_list.append(tile_norm)
        
            if save_tile_file:
                tile_name = "tile_" + str(i_tile).zfill(3) +"_"+ \
                             str(row).zfill(3)+"_" + str(col).zfill(3) + "_" + \
                             "arrayrow_" + str(array_row).zfill(3) + "_col_" + str(array_col).zfill(3) + "_" + \
                             str(downsampling).zfill(3) + "_" + \
                             str(pxl_row_in_fullres) + "_" + str(pxl_col_in_fullres)

                tile_norm.save(f"{tile_folder}/{tile_name}.png")
                tiles_files_list.append(f"{tile_folder}/{tile_name}.png")
        
            ##
            i_tile += 1

    n_tiles = len(tiles_list)
    ##save raw tile data for test
    with open(f'{path2output}/tiles/{slide_name}.norm_tiles_list.pkl', 'wb') as f:
        pickle.dump(tiles_list, f)
else:    
    logger.info(f'Loading file: {path2output}/tiles/{slide_name}.norm_tiles_list.pkl')
    with open(f'{path2output}/tiles/{slide_name}.norm_tiles_list.pkl', 'rb') as f:
        tiles_list = pickle.load(f)
    #with open('tiles/AP67.norm_tiles_list.pkl', 'rb') as f:
    #tiles_list = pickle.load(f)

##====================================================================================================== 
if plot_mask:
    selcted = tissue_potision.iloc[k_barcode_selected]
    selcted.to_csv(f"{path2selected_barcode}/{slide_name}.selected_barcode.csv", index=False)
    ## plot: draw color lines on the mask
    line_color = [20,20,20]
    n_tiles = len(tiles_list)
    img_mask[:,::mask_tile_size,:] = line_color
    img_mask[::mask_tile_size,:,:] = line_color
    mask[:,::mask_tile_size,:] = line_color
    mask[::mask_tile_size,:,:] = line_color
    ##
    fig, ax = plt.subplots(1,2,figsize=(40,20))
    ax[0].imshow(img_mask)
    ax[1].imshow(mask)
    ax[0].set_title(f"{slide_name}, mag_original: {mag_original}, mag_assumed: {mag_assumed}")
    ax[1].set_title(f"n_rows: {n_rows}, n_cols: {n_cols}, n_tiles_total: {n_tiles_total}, n_spot:{tissue_potision.shape[0]+1}, n_tiles_selected: {n_tiles}")
    ##
    plt.tight_layout(h_pad=0.4, w_pad=0.5)
    plt.savefig(f"{path2mask}/{slide_name}_mask.pdf", format="pdf", dpi=50)
    plt.close()
    img_mask = 0 ; mask = 0
    ##
    logger.info("completed cleaning")
    logger.info(f" -- i_slide: {i_slide}, total time: {int(time.time() - start_time)}")

 

n_tiles = len(tiles_list)
logger.info(f"n_tiles: {n_tiles}")
#tiles/AP67.norm_tiles_list.pkl
#with open('tiles/AP67.norm_tiles_list.pkl', 'rb') as f:
#    tiles_list = pickle.load(f)
#n_tiles = len(tiles_list)
##======================================================================================================
##======================================================================================================
batch_size = 64
##----------
## resize:
torch_resize = transforms.Resize(224)

tiles_list_resized = []
for i in range(n_tiles):
    tiles_list_resized.append(torch_resize(tiles_list[i]))

tiles_list_no_resize = tiles_list

tiles_list = tiles_list_resized
tiles_list_resized = []


##----------
## normalize by ImageNet mean and std:
data_transform = transforms.Compose([transforms.ToTensor(),
                                    transforms.Normalize(mean=[0.485, 0.456, 0.406], 
                                                         std=[0.229, 0.224, 0.225])
                                    ])


##----------------------------------------
def extract_features_from_tiles(tiles_list):
    model.eval()
    features = []    
    ## transfor to torch and normalize
    tiles = []
    for i in range(len(tiles_list)):
        tiles.append(data_transform(tiles_list[i]).unsqueeze(0))
    tiles = torch.cat(tiles, dim=0)
    logger.info(f"tiles.shape: {tiles.shape}")   ## [n_tiles, 3, 224, 224]
    #tiles_list = 0
    ##------------------------------------
    ## extract feature from tile image
    features = []
    for idx_start in range(0, n_tiles, batch_size):
        idx_end = idx_start + min(batch_size, n_tiles - idx_start)
        feature = model(tiles[idx_start:idx_end].to(device))        
        features.append(feature.detach().cpu().numpy())
    features = np.concatenate(features)
    return features

##----------------------------------------
if extract_resnet50_features:
    device = (torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu'))
    ## load ResNet pretrained model
    model = Feature_Extraction(model_type="load_from_saved_file")
    model.to(device)
    model.load_state_dict(torch.load(resnet50_file, map_location=device))
    model.eval()
    logger.info("Rsenet50 model loaded")
    features = extract_features_from_tiles(tiles_list)
    logger.info(f"features.shape:, {features.shape}")
    np.save(f"{path2features}{slide_name}_resnet50.npy", features)

## CTransPath
if extract_CTransPath_features:
    logger.info("Extracting CTransPath features\n\n")
    device = torch.device('cpu')
    model = ctranspath()
    model.head = nn.Identity()
    td = torch.load(ctranspath_file)
    model.load_state_dict(td['model'], strict=True)
    model.eval()

    features = extract_features_from_tiles(tiles_list)
    logger.info(f"features.shape:, {features.shape}")
    np.save(f"{path2features}{slide_name}_ctranspath.npy", features)




###MUSK
if extract_MUSK_features:
    ###
    device = (torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu'))
    logger.info("Extracting MUSK features\n\n")    
    tiles_list = tiles_list_no_resize
    n_tiles = len(tiles_list)

    model = create_model("musk_large_patch16_384")
    utils.load_model_and_may_interpolate("hf_hub:xiangjx/musk", model, 'model|module', '')
    model.to(device="cuda", dtype=torch.float16)
    model.eval()

    transform = torchvision.transforms.Compose([
        torchvision.transforms.Resize(384, interpolation=3, antialias=True),
        torchvision.transforms.CenterCrop((384, 384)),
        torchvision.transforms.ToTensor(),
        torchvision.transforms.Normalize(mean=IMAGENET_INCEPTION_MEAN, std=IMAGENET_INCEPTION_STD)
    ])

    def extract_MUSK_features_from_tiles(tiles_list):
        ## transfor to torch and normalize
        features = []
        for i in range(n_tiles):
            img_tensor = transform(tiles_list[i]).unsqueeze(0)        
            with torch.inference_mode():
                image_embeddings = model(
                    image=img_tensor.to(device, dtype=torch.float16),
                    with_head=False,
                    out_norm=False,
                    ms_aug=True,
                    return_global=True  
                )[0]  # return (vision_cls, text_cls)
            features.append(image_embeddings.detach().cpu().numpy()) #2048 features
        features = np.concatenate(features)
        return(features)
    ###
    ###    
    features = extract_MUSK_features_from_tiles(tiles_list) 
    logger.info("features.shape:", features.shape)
    np.save(f"{path2features}{slide_name}_MUSK.npy", features)


###UNI
if extract_UNI_features:
    device = (torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu'))
    logger.info("device:", device)
    logger.info("Extracting UNI features\n\n")
    #with open("hugging_face_token", 'r') as file:
    #    my_token = file.read().strip()
    
    #login(token = my_token)

    local_dir = "/data/wuz6/project/25.pedHGG_spatial/06.Deeplearning/foundataionModels/uni/"
    #os.makedirs(local_dir, exist_ok=True)  # create directory if it does not exist
    #hf_hub_download("MahmoodLab/UNI2-h", filename="pytorch_model.bin", local_dir=local_dir, force_download=True)
    timm_kwargs = {
                'model_name': 'vit_giant_patch14_224',
                'img_size': 224, 
                'patch_size': 14, 
                'depth': 24,
                'num_heads': 24,
                'init_values': 1e-5, 
                'embed_dim': 1536,
                'mlp_ratio': 2.66667*2,
                'num_classes': 0, 
                'no_embed_class': True,
                'mlp_layer': timm.layers.SwiGLUPacked, 
                'act_layer': torch.nn.SiLU, 
                'reg_tokens': 8, 
                'dynamic_img_size': True
            }
    model = timm.create_model(pretrained=False, **timm_kwargs)    
    model.load_state_dict(torch.load(os.path.join(local_dir, "pytorch_model.bin"), map_location= device), strict=True)
    model.eval()

    features = extract_features_from_tiles(tiles_list)
    logger.info("features.shape:", features.shape)
    np.save(f"{path2features}{slide_name}_UNI.npy", features)

    
##======================================================================================================
##======================================================================================================
logger.info(f"finished -- i_slide: {i_slide}, total time: {int(time.time() - start_time)}")





