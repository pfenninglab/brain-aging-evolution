import numpy as np
import pandas as pd
import sys
import os

print(sys.argv[1])
print(sys.argv[2])

PROJ_DIR='/Aging/'+str(sys.argv[1])+'_TACIT/'
CODE_DIR=str(PROJ_DIR)+'code/'
OUT_DIR=str(PROJ_DIR)+'phylolm/'
PREDN_DIR='/machineLearningForComputationalBiology/Cortex_Cell-TACIT/data/tidy_data/240_predictions_matrix_celltypes/'

# Filter OCRs with bh corrected pvalue < 0.05
bh_ocr = pd.read_csv(str(OUT_DIR)+str(sys.argv[2])+'/phylolm_r0_s1_bh_corrected.csv')
bh_ocr = bh_ocr.loc[bh_ocr['OCR'].str.contains('hg38')]
bh_ocr = bh_ocr.loc[bh_ocr['bh']< 0.05]

bh_ocr_pos = bh_ocr.loc[bh_ocr['Coeff']> 0]
bh_ocr_pos = bh_ocr_pos['OCR']
bh_ocr_neg = bh_ocr.loc[bh_ocr['Coeff']< 0]
bh_ocr_neg = bh_ocr_neg['OCR']

# Split OCR into chr, start and end positions
# Function to process each string
def process_string(data):
    split_data = data.split(':')
    chr = split_data[1]
    positions = split_data[2].split('-')
    start = int(positions[0])
    end = int(positions[1])
    return [chr, start, end]

# Apply the function to each string in the list
pos = [process_string(data) for data in bh_ocr_pos]
pos = pd.DataFrame(pos, columns=['chr', 'start', 'end'])
neg = [process_string(data) for data in bh_ocr_neg]
neg = pd.DataFrame(neg, columns=['chr', 'start', 'end'])

# Save as bed files
pos.to_csv(str(OUT_DIR)+str(sys.argv[2])+'/phylolm_r0_s1_bh_corrected_pos.bed', index = False, sep ='\t', header = None)
neg.to_csv(str(OUT_DIR)+str(sys.argv[2])+'/phylolm_r0_s1_bh_corrected_neg.bed', index = False, sep ='\t', header = None)

