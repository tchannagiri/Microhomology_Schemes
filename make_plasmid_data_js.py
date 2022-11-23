import json
import pandas as pd
import os

MMEJ_LIST_CSV = os.path.join('csv', 'plasmid', 'Sense_MH_list.csv')
MMEJ_LIST_ANTI_CSV = os.path.join('csv', 'plasmid', 'Antisense_MH_list.csv')
SCHEMES_DATA_JS = os.path.join('html', 'plasmid', 'data.js')
# MICROHOMOLOGIES_VAR_JS = 'MICROHOMOLOGIES'
# AREAS_VAR_JS = 'AREAS'
# BARS_VAR_JS = 'BARS'
# PRIMER_VAR_JS = 'PRIMER'
# CUT_POS_VAR_JS = 'CUT_POS'
# REF_SEQ_VAR_JS = 'REF_SEQ'

JS_VARS = {
  'REF_SEQ': {
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
  },
  'AREAS': {
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
    '2dsb': {
        'Forward': [ 
          {'name': 'exon1', 'start': 21, 'end': 66},
          {'name': 'intron', 'start': 67, 'end': 195},
          {'name': 'splice', 'start': 196, 'end': 201},
          {'name': 'exon2', 'start': 202, 'end': 233},
        ],
        'Reverse': [
          {'name': 'exon2', 'start': 21, 'end': 52},
          {'name': 'splice', 'start': 53, 'end': 58},
          {'name': 'intron', 'start': 59, 'end': 187},
          {'name': 'exon1', 'start': 188, 'end': 233},
        ],
      },
    },
    'd5': {
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
  },
  'PRIMER': {
    'wt': {
      'Forward': [
        {'name': 'p1', 'start': 1, 'end': 20},
        {'name': 'p2', 'start': 210, 'end': 229},
      ],
      'Reverse': [
        {'name': 'p2', 'start': 1, 'end': 20},
        {'name': 'p1', 'start': 210, 'end': 229},
      ],
    },
    'db': {
      'Forward': [
        {'name': 'p1', 'start': 1, 'end': 20},
        {'name': 'p2', 'start': 155, 'end': 174},
      ],
      'Reverse': [
        {'name': 'p2', 'start': 1, 'end': 20},
        {'name': 'p1', 'start': 155, 'end': 174},
      ],
    },
    'awt': {
      'Forward': [
        {'name': 'p1', 'start': 1, 'end': 20},
        {'name': 'p2', 'start': 234, 'end': 253},
      ],
      'Reverse': [
        {'name': 'p2', 'start': 1, 'end': 20},
        {'name': 'p1', 'start': 234, 'end': 253},
      ],
    },
    'd5': {
      'Forward': [
        {'name': 'p1', 'start': 1, 'end': 20},
        {'name': 'p2', 'start': 228, 'end': 247},
      ],
      'Reverse': [
        {'name': 'p2', 'start': 1, 'end': 20},
        {'name': 'p1', 'start': 228, 'end': 247},
      ],
    },
  },
  'BARS': {
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
  },
  'CUT_POS': {
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
  },
  'BAR_COLOR': {
    'dsred': '#4472c4',
    'intron': '#00b050',
  },
  'BAR_TEXT': {
    'dsred': 'DS Red',
    'intron': 'Intron',
  },
  'AREAS_COLOR': {
    'exon1': '#4472c4', 
    'exon2': '#4472c4',
    'intron': '#00b050',
    'branch': '#548235',
    'splice': '#548235', 
  },
  'AREAS_TEXT': {
    'exon1': 'Exon1',
    'exon2': 'Exon2',
    'intron': 'Intron',
    'branch': 'Intron (branch site)',
    'splice': '5\' splice',
  },
  'PRIMER_TEXT': {
    'p2': 'Primer',
    'p1': 'Primer',
  },
  'LABELS': {
    'hg39': 'sgRNA A',
    'hg42': 'sgRNA B',
    'Forward': 'Forward strand',
    'Reverse': 'Reverse strand',
    'wt': 'Sense/pCMVΔ',
    'db': 'BranchΔ',
    'awt': 'Antisense',
    'd5': '5\' splicingΔ',
    'all': None,
    'exon_exon': 'Exon-Exon group',
    'exon_intron': 'Exon-Intron group',
  },
  'CELLTYPE_FLIPPED_INTRON': ['d5', 'awt'],
}

microhomologies = [
  pd.read_csv(MMEJ_LIST_CSV),
  pd.read_csv(MMEJ_LIST_ANTI_CSV),
]
microhomologies = pd.concat(microhomologies, axis='index')
microhomologies = microhomologies[['Name', 'Celltype', 'Breaks', 'Type', 'Strand', 'Left', 'Right', 'Pattern', 'Bold']]
microhomologies = microhomologies.to_dict('records')

with open(SCHEMES_DATA_JS, 'w') as out:
  out.write(f'var MICROHOMOLOGIES = ')
  json.dump(microhomologies, out, indent=2)
  out.write(';\n\n')

  for name, value in JS_VARS.items():
    out.write(f'var {name} = ')
    json.dump(value, out, indent=2)
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

