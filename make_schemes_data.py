import json
import pandas as pd

MMEJ_LIST_CSV = 'Sense_MH_list.csv'
MMEJ_LIST_ANTI_CSV = 'Antisense_MH_list.csv'
SCHEMES_DATA_JS = 'schemes_data.js'
MICROHOMOLOGIES_VAR_JS = 'MICROHOMOLOGIES'
AREAS_VAR_JS = 'AREAS'
BARS_VAR_JS = 'BARS'
PCIS_VAR_JS = 'PCIS'
CUT_POS_VAR_JS = 'CUT_POS'
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
  'awt': {
    'Forward': 'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTTAGGGATAACAGGGTAATACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
    'Reverse': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTATTACCCTGTTATCCCTAACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
  },
  'd5': {
    'Forward': 'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTTAGGGATAACAGGGTAATACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTACTGCAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
    'Reverse': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGCAGTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTATTACCCTGTTATCCCTAACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
  },
  # For the OLD Antisense reference sequences
  # 'awt': {
  #   'Forward': 'TGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTTAGGGATAACAGGGTAATACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
  #   'Reverse': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTATTACCCTGTTATCCCTAACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCA',
  # },
  # 'd5': {
  #   'Forward': 'TGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTTAGGGATAACAGGGTAATACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTACTGCAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
  #   'Reverse': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGCAGTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTATTACCCTGTTATCCCTAACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCA',
  # },
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
        {'name': 'exon1', 'start': 158, 'end': 209},
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
        {'name': 'exon1', 'start': 103, 'end': 154},
      ],
    },
  },
  'awt': {
   '2dsb': { # len 253
      'Forward': [ 
        {'name': 'exon1', 'start': 21, 'end': 66},
        {'name': 'intron', 'start': 67, 'end': 195},
        {'name': 'splice', 'start': 196, 'end': 201},
        {'name': 'exon2', 'start': 202, 'end': 233},
      ],
      'Reverse': [ # len 253
        {'name': 'exon2', 'start': 21, 'end': 52},
        {'name': 'splice', 'start': 53, 'end': 58},
        {'name': 'intron', 'start': 59, 'end': 187},
        {'name': 'exon1', 'start': 188, 'end': 233},
      ],
    },
  },
  'd5': { # len 247
   '2dsb': {
      'Forward': [
        {'name': 'exon1', 'start': 21, 'end': 66},
        {'name': 'intron', 'start': 67, 'end': 195},
        {'name': 'exon2', 'start': 196, 'end': 227},
      ],
      'Reverse': [
        {'name': 'exon2', 'start': 21, 'end': 52},
        {'name': 'intron', 'start': 53, 'end': 181},
        {'name': 'exon1', 'start': 182, 'end': 227},
      ],
    },
  },
  # For the OLD Antisense reference sequences
  # 'awt': {
  #  '2dsb': {
  #     'Forward': [ #len 237
  #       {'name': 'exon1', 'start': 21, 'end': 50},
  #       {'name': 'intron', 'start': 51, 'end': 179},
  #       {'name': 'splice', 'start': 180, 'end': 185},
  #       {'name': 'exon2', 'start': 186, 'end': 217},
  #     ],
  #     'Reverse': [
  #       {'name': 'exon2', 'start': 21, 'end': 52},
  #       {'name': 'splice', 'start': 53, 'end': 58},
  #       {'name': 'intron', 'start': 59, 'end': 187},
  #       {'name': 'exon1', 'start': 188, 'end': 217},
  #     ],
  #   },
  # },
  # 'd5': { # len 231
  #  '2dsb': {
  #     'Forward': [
  #       {'name': 'exon1', 'start': 21, 'end': 50},
  #       {'name': 'intron', 'start': 51, 'end': 179},
  #       {'name': 'exon2', 'start': 180, 'end': 211},
  #     ],
  #     'Reverse': [
  #       {'name': 'exon2', 'start': 21, 'end': 52},
  #       {'name': 'intron', 'start': 53, 'end': 181},
  #       {'name': 'exon1', 'start': 182, 'end': 211},
  #     ],
  #   },
  # },
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
  },
  'awt': {
    'Forward': [
      {'name': 'f6', 'start': 1, 'end': 20},
      {'name': 'r5', 'start': 234, 'end': 253},
    ],
    'Reverse': [
      {'name': 'r5', 'start': 1, 'end': 20},
      {'name': 'f6', 'start': 234, 'end': 253},
    ],
  },
  'd5': {
    'Forward': [
      {'name': 'f6', 'start': 1, 'end': 20},
      {'name': 'r5', 'start': 228, 'end': 247},
    ],
    'Reverse': [
      {'name': 'r5', 'start': 1, 'end': 20},
      {'name': 'f6', 'start': 228, 'end': 247},
    ],
  },
  # For the OLD Antisense reference sequences
  # 'awt': {
  #   'Forward': [
  #     {'name': 'f6', 'start': 1, 'end': 20},
  #     {'name': 'r5', 'start': 218, 'end': 237},
  #   ],
  #   'Reverse': [
  #     {'name': 'r5', 'start': 1, 'end': 20},
  #     {'name': 'f6', 'start': 218, 'end': 237},
  #   ],
  # },
  # 'd5': {
  #   'Forward': [
  #     {'name': 'f6', 'start': 1, 'end': 20},
  #     {'name': 'r5', 'start': 212, 'end': 231},
  #   ],
  #   'Reverse': [
  #     {'name': 'r5', 'start': 1, 'end': 20},
  #     {'name': 'f6', 'start': 212, 'end': 231},
  #   ],
  # },
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
  'awt': {
    'Forward': [
      {'name': 'dsred', 'start': 1, 'end': 66},
      {'name': 'intron', 'start': 67, 'end': 201},
      {'name': 'dsred', 'start': 202, 'end': 253},
    ],
    'Reverse': [
      {'name': 'dsred', 'start': 1, 'end': 52},
      {'name': 'intron', 'start': 53, 'end': 187},
      {'name': 'dsred', 'start': 188, 'end': 253},
    ],
  },
  'd5': {
    'Forward': [
      {'name': 'dsred', 'start': 1, 'end': 66},
      {'name': 'intron', 'start': 67, 'end': 195},
      {'name': 'dsred', 'start': 196, 'end': 247},
    ],
    'Reverse': [
      {'name': 'dsred', 'start': 1, 'end': 52},
      {'name': 'intron', 'start': 53, 'end': 181},
      {'name': 'dsred', 'start': 182, 'end': 247},
    ],
  },
  # For the OLD Antisense reference sequences
  # 'awt': {
  #   'Forward': [
  #     {'name': 'dsred', 'start': 1, 'end': 50},
  #     {'name': 'intron', 'start': 51, 'end': 185},
  #     {'name': 'dsred', 'start': 186, 'end': 237},
  #   ],
  #   'Reverse': [
  #     {'name': 'dsred', 'start': 1, 'end': 52},
  #     {'name': 'intron', 'start': 53, 'end': 187},
  #     {'name': 'dsred', 'start': 188, 'end': 237},
  #   ],
  # },
  # 'd5': {
  #   'Forward': [
  #     {'name': 'dsred', 'start': 1, 'end': 50},
  #     {'name': 'intron', 'start': 51, 'end': 179},
  #     {'name': 'dsred', 'start': 180, 'end': 231},
  #   ],
  #   'Reverse': [
  #     {'name': 'dsred', 'start': 1, 'end': 52},
  #     {'name': 'intron', 'start': 53, 'end': 181},
  #     {'name': 'dsred', 'start': 182, 'end': 231},
  #   ],
  # },
}

CUT_POS = {
  'wt': {
    'hg39': {
      'Forward': [67], 
      'Reverse':  [162],
    },
    'hg42': {
      'Forward': [183], 
      'Reverse':  [46],
    },
    '2dsb': {
      'Forward': [67, 183],
      'Reverse': [46, 162],
    },
  },
  'db': {
    'hg39': {
      'Forward': [67], 
      'Reverse':  [107],
    },
    'hg42': {
      'Forward': [128], 
      'Reverse':  [46],
    },
    '2dsb': {
      'Forward': [67, 128],
      'Reverse': [46, 107],
    },
  },
  'awt': {
    '2dsb': {
      'Forward': [66, 206],
      'Reverse': [47, 187],
    },
  },
  'd5': {
    '2dsb': {
      'Forward': [66, 200],
      'Reverse': [47, 181],
    },
  },
  # For the OLD Antisense reference sequences
  # 'awt': {
  #   '2dsb': {
  #     'Forward': [50, 190],
  #     'Reverse': [47, 187],
  #   },
  # },
  # 'd5': {
  #   '2dsb': {
  #     'Forward': [50, 184],
  #     'Reverse': [47, 181],
  #   },
  # },
}

microhomologies = [
  pd.read_csv(MMEJ_LIST_CSV),
  pd.read_csv(MMEJ_LIST_ANTI_CSV),
]
microhomologies = pd.concat(microhomologies, axis='index')
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
  
  out.write(f'var {CUT_POS_VAR_JS} = ')
  json.dump(CUT_POS, out, indent=2)
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

