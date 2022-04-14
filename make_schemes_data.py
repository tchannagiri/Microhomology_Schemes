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

MMEJ_LIST_CSV = 'MMEJ_list_bold.csv'
SCHEMES_DATA_JS = 'schemes_data.js'
MICROHOMOLOGIES_VAR_JS = 'MICROHOMOLOGIES'
AREAS_VAR_JS = 'AREAS'
BARS_VAR_JS = 'BARS'
PCIS_VAR_JS = 'PCIS'
REF_SEQ_VAR_JS = 'REF_SEQ'

REF_SEQ = {
  'wt': {
    'Forward': 'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
    'Reverse': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
  },
  'db': {
    'Forward': 'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
    'Reverse': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGTGGATAACTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
  },
}

AREAS = {
  'wt': {
    'hg39': {
      'Forward': [
        {'name': 'exon1', 'start': 21, 'end': 72},
        {'name': 'intron', 'start': 73, 'end': 118},
        {'name': 'branch', 'start': 119, 'end': 173},
        {'name': 'intron', 'start': 174, 'end': 183},
        {'name': 'exon2', 'start': 184, 'end': 209},
      ],
      'Reverse': [
        {'name': 'exon2', 'start': 21, 'end': 46},
        {'name': 'intron', 'start': 47, 'end': 56},
        {'name': 'branch', 'start': 57, 'end': 111},
        {'name': 'intron', 'start': 112, 'end': 157},
        {'name': 'exon1', 'start': 158, 'end': 209},
      ],
    },
    'hg42': {
      'Forward': [
        {'name': 'exon1', 'start': 21, 'end': 72},
        {'name': 'intron', 'start': 73, 'end': 118},
        {'name': 'branch', 'start': 119, 'end': 173},
        {'name': 'intron', 'start': 174, 'end': 183},
        {'name': 'exon2', 'start': 184, 'end': 209},
      ],
      'Reverse': [
        {'name': 'exon2', 'start': 21, 'end': 46},
        {'name': 'intron', 'start': 47, 'end': 56},
        {'name': 'branch', 'start': 57, 'end': 111},
        {'name': 'intron', 'start': 112, 'end': 157},
        {'name': 'exon1', 'start': 158, 'end': 209},
      ],
    },
    '2dsb': {
      'Forward': [
        {'name': 'exon1', 'start': 21, 'end': 72},
        {'name': 'intron', 'start': 73, 'end': 118},
        {'name': 'branch', 'start': 119, 'end': 173},
        {'name': 'intron', 'start': 174, 'end': 183},
        {'name': 'exon2', 'start': 184, 'end': 209},
      ],
      'Reverse': [
        {'name': 'exon2', 'start': 21, 'end': 46},
        {'name': 'intron', 'start': 47, 'end': 56},
        {'name': 'branch', 'start': 57, 'end': 111},
        {'name': 'intron', 'start': 112, 'end': 157},
        {'name': 'exon1', 'start': 163, 'end': 209},
      ],
    },
  },
  'db': {
   'hg39': {
      'Forward': [
        {'name': 'exon1', 'start': 21, 'end': 72},
        {'name': 'intron', 'start': 73, 'end': 128},
        {'name': 'exon2', 'start': 129, 'end': 154},
      ],
      'Reverse': [
        {'name': 'exon2', 'start': 21, 'end': 46},
        {'name': 'intron', 'start': 47, 'end': 102},
        {'name': 'exon1', 'start': 103, 'end': 154},
      ],
    },
    'hg42': {
      'Forward': [
        {'name': 'exon1', 'start': 21, 'end': 72},
        {'name': 'intron', 'start': 73, 'end': 128},
        {'name': 'exon2', 'start': 129, 'end': 154},
      ],
      'Reverse': [
        {'name': 'exon2', 'start': 21, 'end': 46},
        {'name': 'intron', 'start': 47, 'end': 102},
        {'name': 'exon1', 'start': 103, 'end': 154},
      ],
    },
    '2dsb': {
      'Forward': [
        {'name': 'exon1', 'start': 21, 'end': 72},
        {'name': 'intron', 'start': 73, 'end': 128},
        {'name': 'exon2', 'start': 129, 'end': 154},
      ],
      'Reverse': [
        {'name': 'exon2', 'start': 21, 'end': 46},
        {'name': 'intron', 'start': 47, 'end': 102},
        {'name': 'exon1', 'start': 108, 'end': 154},
      ],
    },
  },
}

PCIS = {
  'wt': {
    'Forward': [
      {'name': 'f6', 'start': 1, 'end': 20},
      {'name': 'r5', 'start': 210, 'end': 229},
    ],
    'Reverse': [
      {'name': 'r5', 'start': 1, 'end': 20},
      {'name': 'f6', 'start': 210, 'end': 229},
    ],
  },
  'db': {
    'Forward': [
      {'name': 'f6', 'start': 1, 'end': 20},
      {'name': 'r5', 'start': 155, 'end': 174},
    ],
    'Reverse': [
      {'name': 'r5', 'start': 1, 'end': 20},
      {'name': 'f6', 'start': 155, 'end': 174},
    ],
  }
}

BARS = {
  'wt': {
    'Forward': [
      {'name': 'dsred', 'start': 1, 'end': 72},
      {'name': 'intron', 'start': 73, 'end': 183},
      {'name': 'dsred', 'start': 184, 'end': 229},
    ],
    'Reverse': [
      {'name': 'dsred', 'start': 1, 'end': 46},
      {'name': 'intron', 'start': 47, 'end': 157},
      {'name': 'dsred', 'start': 158, 'end': 229},
    ],
  },
  'db': {
    'Forward': [
      {'name': 'dsred', 'start': 1, 'end': 72},
      {'name': 'intron', 'start': 73, 'end': 128},
      {'name': 'dsred', 'start': 129, 'end': 174},
    ],
    'Reverse': [
      {'name': 'dsred', 'start': 1, 'end': 46},
      {'name': 'intron', 'start': 47, 'end': 102},
      {'name': 'dsred', 'start': 103, 'end': 174},
    ],
  },
}

microhomologies = pd.read_csv(MMEJ_LIST_CSV)
microhomologies = microhomologies[['Name', 'Celltype', 'Breaks', 'Type', 'Strand', 'Left', 'Right', 'Pattern', 'Bold']]
microhomologies = microhomologies.to_dict('records')

with open(SCHEMES_DATA_JS, 'w') as out:
  out.write(f'var {MICROHOMOLOGIES_VAR_JS} = ')
  json.dump(microhomologies, out, indent=2)
  out.write(';\n\n')

  out.write(f'var {AREAS_VAR_JS} = ')
  json.dump(AREAS, out, indent=2)
  out.write(';\n\n')
  
  out.write(f'var {BARS_VAR_JS} = ')
  json.dump(BARS, out, indent=2)
  out.write(';\n\n')
  
  out.write(f'var {PCIS_VAR_JS} = ')
  json.dump(PCIS, out, indent=2)
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


