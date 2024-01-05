import json
import pandas as pd
import os

MMEJ_LIST_CSV = [
  os.path.join('csv', 'plasmid_2', 'mmej_color_sgRNAE_1DSB.csv'),
  os.path.join('csv', 'plasmid_2', 'mmej_color_sgRNAE_sgRNAI_2DSB.csv'),
  os.path.join('csv', 'plasmid_2', 'mmej_color_sgRNAE_sgRNAJ_2DSB.csv'),
  os.path.join('csv', 'plasmid_2', 'mmej_color_sgRNAK_sgRNAL_2DSB_Exon_Exon.csv'),
  os.path.join('csv', 'plasmid_2', 'mmej_color_sgRNAK_sgRNAL_2DSB_intron_intron.csv'),
]
SCHEMES_DATA_JS = os.path.join('html', 'plasmid_2', 'data.js')

JS_VARS = {
  'BOLDED': {
    'Default': {
      'wt': {
        'E': [
          'EE1', 'EE3,' 'EE9,' 'EE10', 'EE11', 'EE12', 'EE15', 'EE17', 'EE20', 'EE24', 'EE28',
          'EE35', 'EE36', 'EE38', 'EE40', 'EE41', 'EE43', 'EE44', 'EE45', 'EE46', 'EE47', 'EE51',
          'EE53', 'EE54', 'EE55', 'EE57', 'EE58', 'EE61', 'EE64', 'EE65', 'EE67', 'EE68', 'EE71',
          'EE73', 'EE74', 'EE75', 'EE76', 'EE77', 'EE78', 'EE79', 'EE80', 'EE81', 'EE82', 'EE85',
          'EE86', 'EE87', 'EE88', 'EE89', 'EE90', 'EE91', 'EE92', 'EE93', 'EE94', 'EE95', 'EE96',
          'EE97', 'EE98', 'EE99', 'EE100', 'EE101',
        ],
        'EI': [],
        'EJ': [
          'EE1', 'EE2', 'EE4', 'EE5', 'EE6', 'EE7', 'EE8', 'EE9', 'EE10', 'EE11', 'EE12', 'EE13',
          'EE14', 'EE15', 'EE16', 'EE17', 'EE18', 'EE19', 'EE20', 'EE21', 'EE22', 'EE23',
          'EE1R', 'EE2R', 'EE4R', 'EE5R', 'EE6R', 'EE7R', 'EE8R', 'EE9R', 'EE10R', 'EE11R', 'EE12R', 'EE13R',
          'EE14R', 'EE15R', 'EE16R', 'EE17R', 'EE18R', 'EE19R', 'EE20R', 'EE21R', 'EE22R', 'EE23R',
        ],
      },
      'db': {
        'E': [
          'EE1', 'EE3,' 'EE9,' 'EE10', 'EE11', 'EE12', 'EE15', 'EE17', 'EE20', 'EE24', 'EE28',
          'EE35', 'EE36', 'EE38', 'EE40', 'EE41', 'EE43', 'EE44', 'EE45', 'EE46', 'EE47', 'EE51',
          'EE53', 'EE54', 'EE55', 'EE57', 'EE58', 'EE61', 'EE64', 'EE65', 'EE67', 'EE68', 'EE71',
          'EE73', 'EE74', 'EE75', 'EE76', 'EE77', 'EE78', 'EE79', 'EE80', 'EE81', 'EE82', 'EE85',
          'EE86', 'EE87', 'EE88', 'EE89', 'EE90', 'EE91', 'EE92', 'EE93', 'EE94', 'EE95', 'EE96',
          'EE97', 'EE98', 'EE99', 'EE100', 'EE101',
        ],
        'EI': [],
        'EJ': [
          'EE1', 'EE2', 'EE4', 'EE5', 'EE6', 'EE7', 'EE8', 'EE9', 'EE10', 'EE11', 'EE12', 'EE13',
          'EE14', 'EE15', 'EE16', 'EE17', 'EE18', 'EE19', 'EE20', 'EE21', 'EE22', 'EE23',
          'EE1R', 'EE2R', 'EE4R', 'EE5R', 'EE6R', 'EE7R', 'EE8R', 'EE9R', 'EE10R', 'EE11R', 'EE12R', 'EE13R',
          'EE14R', 'EE15R', 'EE16R', 'EE17R', 'EE18R', 'EE19R', 'EE20R', 'EE21R', 'EE22R', 'EE23R',
        ],
      },
      'awt': {
        'KL': [],
      },
      'd5': {
        'KL': [],
      },
    }
  },
  'REF_SEQ': {
    'wt': {
      'Forward': 'CTTCAAGGTGCGCATGGAGGGCACCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCCACAACACCGTGAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCCCAGTTCCAGTACGGCTCCAAGGTGTACGTGAAGCACCCCGCCGACATCCCCGACTACAAGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGAT',
      'Reverse': 'ATCACGCGCTCCCACTTGAAGCCCTCGGGGAAGGACAGCTTCTTGTAGTCGGGGATGTCGGCGGGGTGCTTCACGTACACCTTGGAGCCGTACTGGAACTGGGGGGACAGGATGTCCCAGGCGAAGGGCAGGGGGCCGCCCTTGGTCACCTTCAGCTTCACGGTGTTGTGGCCCTCGTAGGGGCGGCCCTCGCCCTCGCCCTCGATCTCGAACTCGTGGCCGTTCACGGTGCCCTCCATGCGCACCTTGAAG',
    },
    'db': {
      'Forward': 'CTTCAAGGTGCGCATGGAGGGCACCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCCACAACACCGTGAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCCCAGTTCCAGTACGGCTCCAAGGTGTACGTGAAGCACCCCGCCGACATCCCCGACTACAAGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGAT',
      'Reverse': 'ATCACGCGCTCCCACTTGAAGCCCTCGGGGAAGGACAGCTTCTTGTAGTCGGGGATGTCGGCGGGGTGCTTCACGTACACCTTGGAGCCGTACTGGAACTGGGGGGACAGGATGTCCCAGGCGAAGGGCAGGGGGCCGCCCTTGGTCACCTTCAGCTTCACGGTGTTGTGGCCCTCGTAGGGGCGGCCCTCGCCCTCGCCCTCGATCTCGAACTCGTGGCCGTTCACGGTGCCCTCCATGCGCACCTTGAAG',
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
      'E': {
        'Forward': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
        'Reverse': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
      },
      'EI': {
        'Forward': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
        'Reverse': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
      },
      'EJ': {
        'Forward': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
        'Reverse': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
      },
    },
    'db': {
      'E': {
        'Forward': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
        'Reverse': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
      },
      'EI': {
        'Forward': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
        'Reverse': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
      },
      'EJ': {
        'Forward': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
        'Reverse': [
          {'name': 'exon1', 'start': 21, 'end': 232},
        ],
      },
    },
    'awt': {
      'KL': {
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
      }
    },
    'd5': {
      'KL': {
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
      }
    },
  },
  'PRIMER': {
    'wt': {
      'Forward': [
        {'name': 'p1', 'start': 1, 'end': 20},
        {'name': 'p2', 'start': 233, 'end': 252},
      ],
      'Reverse': [
        {'name': 'p2', 'start': 1, 'end': 20},
        {'name': 'p1', 'start': 233, 'end': 252},
      ],
    },
    'db': {
      'Forward': [
        {'name': 'p1', 'start': 1, 'end': 20},
        {'name': 'p2', 'start': 233, 'end': 252},
      ],
      'Reverse': [
        {'name': 'p2', 'start': 1, 'end': 20},
        {'name': 'p1', 'start': 233, 'end': 252},
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
        {'name': 'dsred', 'start': 1, 'end': 252},
      ],
      'Reverse': [
        {'name': 'dsred', 'start': 1, 'end': 252},
      ],
    },
    'db': {
      'Forward': [
        {'name': 'dsred', 'start': 1, 'end': 252},
      ],
      'Reverse': [
        {'name': 'dsred', 'start': 1, 'end': 252},
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
      'E': {
        'Forward': [72], 
        'Reverse':  [180],
      },
      'EI': {
        'Forward': [72, 130], 
        'Reverse':  [122, 180],
      },
      'EJ': {
        'Forward': [72, 194], 
        'Reverse':  [58, 180],
      },
    },
    'db': {
      'E': {
        'Forward': [72], 
        'Reverse':  [180],
      },
      'EI': {
        'Forward': [72, 130], 
        'Reverse':  [122, 180],
      },
      'EJ': {
        'Forward': [72, 194], 
        'Reverse':  [58, 180],
      },
    },
    'awt': {
      'KL': {
        'Forward': [111, 138], 
        'Reverse':  [115, 142],
      },
    },
    'd5': {
      'KL': {
        'Forward': [111, 138], 
        'Reverse':  [109, 136],
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
    'exon1': '#ff0000',
    'exon2': '#ff0000',
    'intron': '#548235',
    'splice': '#548235',
  },
  'AREAS_TEXT': {
    'exon1': 'Exon1',
    'exon2': 'Exon2',
    'intron': 'Intron',
    'splice': '5\' splice',
  },
  'PRIMER_TEXT': {
    'p1': 'Primer',
    'p2': 'Primer',
  },
  'LABELS': {
    'E': 'sgRNA E',
    'EI': 'sgRNA E & I',
    'EJ': 'sgRNA E & J',
    'KL': 'sgRNA E & K',
    'Forward': 'forward strand',
    'Reverse': 'reverse strand',
    'wt': 'Sense/pCMVΔ',
    'db': 'BranchΔ',
    'awt': 'Antisense',
    'd5': '5\' SplicingΔ',
    'all': None,
    'exon_exon': 'Exon-Exon group',
    'intron_intron': 'Intron-Intron group',
  },
  'CELLTYPE_FLIPPED_INTRON': ['d5', 'awt'],
}

microhomologies = [pd.read_csv(x) for x in MMEJ_LIST_CSV]
microhomologies = pd.concat(microhomologies, axis='index')
microhomologies = microhomologies[['Name', 'Celltype', 'Breaks', 'Type', 'Strand', 'Left', 'Right', 'Pattern']]
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
#     "microhomology_id": "spliced_HG39_p2_G1_exon1_exon2_P21_P207_GAA",
#     "hguide": "HG39",
#     "strand": "p2",
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

