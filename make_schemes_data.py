import json
import time
import os
import pathlib
import pandas as pd
import numpy as np
import re
import seaborn as sns
import random
import itertools

MMEJ_LIST_CSV = 'data/MMEJ_list.csv'
SCHEMES_DATA_JS = 'schemes_data.js'
MICROHOMOLOGIES_VAR_JS = 'MICROHOMOLOGIES'
AREAS_VAR_JS = 'AREAS'
REF_SEQ_VAR_JS = 'REF_SEQ'

REF_SEQ = {
  'Forward': 'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
  'Reverse': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
}

AREAS = {
  'hg39': {
    'Forward': [
      {'name': 'exon1', 'start': 21, 'end': 67},
      {'name': 'exon2', 'start': 184, 'end': 209},
      {'name': 'intron', 'start': 73, 'end': 118},
      {'name': 'branch', 'start': 119, 'end': 173},
      {'name': 'intron', 'start': 174, 'end': 183},
    ],
    'Reverse': [
      {'name': 'exon1', 'start': 163, 'end': 209},
      {'name': 'exon2', 'start': 21, 'end': 46},
      {'name': 'intron', 'start': 47, 'end': 56},
      {'name': 'branch', 'start': 57, 'end': 111},
      {'name': 'intron', 'start': 112, 'end': 157},
    ],
  },
  'hg42': {
    'Forward': [
      {'name': 'exon2', 'start': 21, 'end': 46},
      {'name': 'exon1', 'start': 158, 'end': 209},
      {'name': 'intron', 'start': 47, 'end': 56},
      {'name': 'branch', 'start': 57, 'end': 111},
      {'name': 'intron', 'start': 112, 'end': 157},
    ],
    'Reverse': [
      {'name': 'exon2', 'start': 184, 'end': 209},
      {'name': 'exon1', 'start': 21, 'end': 72},
      {'name': 'intron', 'start': 73, 'end': 118},
      {'name': 'branch', 'start': 119, 'end': 173},
      {'name': 'intron', 'start': 174, 'end': 183},
    ],
  }
}

microhomologies = pd.read_csv(MMEJ_LIST_CSV)
microhomologies = microhomologies[['Name', 'Celltype', 'Breaks', 'Type', 'Strand', 'Left', 'Right', 'Pattern']]
microhomologies = microhomologies.to_dict('records')

with open(SCHEMES_DATA_JS, 'w') as out:
  out.write(f'var {MICROHOMOLOGIES_VAR_JS} = ')
  json.dump(microhomologies, out, indent=2)
  out.write(';\n\n')

  out.write(f'var {AREAS_VAR_JS} = ')
  json.dump(AREAS, out, indent=2)
  out.write(';\n\n')
  
  out.write(f'var {REF_SEQ_VAR_JS} = ')
  json.dump(REF_SEQ, out, indent=2)
  out.write(';\n\n')



#  {
#     "microhomology_id": "spliced_HG39_R1_G1_exon1_exon2_P21_P207_GAA",
#     "hguide": "HG39",
#     "strand": "R1",
#     "treatment": "spliced",
#     "microhomology_group": "G1",
#     "match_length": 3,
#     "microhomology_match": "GAA",
#     "pos_left": 21,
#     "pos_right": 207,
#     "microhomology_repaired": "TTCAAGTGGGAGCGCGTGATGAAGTTCATCGGCGTGAACTTCC",
#     "window_around_repair": "AGCGCGTGATGAAGTTCATCGGC",
#     "window_radius": 10,
#     "ref_seq": "TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC",
#     "visual_comparison": "++++++++++++++++++++MMM---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------MMM++++++++++++++++++++",
#     "scheme": "                    EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE                                                                                                                    EEEEEEEEEEEEEEEEEEEEEEEEEE                    "
#   },


