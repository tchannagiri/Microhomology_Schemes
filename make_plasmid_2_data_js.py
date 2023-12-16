import json
import pandas as pd
import os

MMEJ_LIST_CSV = [
  os.path.join('csv', 'plasmid_2', 'mmej_color_sgRNAE_1DSB.csv'),
  os.path.join('csv', 'plasmid_2', 'mmej_color_sgRNAE_sgRNAI_2DSB.csv'),
  os.path.join('csv', 'plasmid_2', 'mmej_color_sgRNAE_sgRNAJ_2DSB.csv'),
  # MAKE THE CSV FILES FOR Antisense and 5'-Splicing!!!
]
SCHEMES_DATA_JS = os.path.join('html', 'plasmid_2', 'data.js')

JS_VARS = {
  'BOLDED': {
    'Default': {
      'wt': {
        'E': [],
        'EI': [],
        'EJ': [],
      },
      'db': {
        'E': [],
        'EI': [],
        'EJ': [],
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
          {'name': 'region1', 'start': 21, 'end': 72},
          {'name': 'region2', 'start': 73, 'end': 232},
        ],
        'Reverse': [
          {'name': 'region2', 'start': 21, 'end': 180},
          {'name': 'region1', 'start': 181, 'end': 232},
        ],
      },
      'EI': {
        'Forward': [
          {'name': 'region1', 'start': 21, 'end': 72},
          {'name': 'region2', 'start': 131, 'end': 232},
        ],
        'Reverse': [
          {'name': 'region2', 'start': 21, 'end': 122},
          {'name': 'region1', 'start': 181, 'end': 232},
        ],
      },
      'EJ': {
        'Forward': [
          {'name': 'region1', 'start': 21, 'end': 72},
          {'name': 'region2', 'start': 195, 'end': 232},
        ],
        'Reverse': [
          {'name': 'region2', 'start': 21, 'end': 58},
          {'name': 'region1', 'start': 181, 'end': 232},
        ],
      },
    },
    'db': {
      'E': {
        'Forward': [
          {'name': 'region1', 'start': 21, 'end': 72},
          {'name': 'region2', 'start': 73, 'end': 232},
        ],
        'Reverse': [
          {'name': 'region2', 'start': 1, 'end': 160},
          {'name': 'region1', 'start': 161, 'end': 232},
        ],
      },
      'EI': {
        'Forward': [
          {'name': 'region1', 'start': 21, 'end': 72},
          {'name': 'region2', 'start': 131, 'end': 232},
        ],
        'Reverse': [
          {'name': 'region2', 'start': 21, 'end': 122},
          {'name': 'region1', 'start': 181, 'end': 232},
        ],
      },
      'EJ': {
        'Forward': [
          {'name': 'region1', 'start': 21, 'end': 72},
          {'name': 'region2', 'start': 195, 'end': 232},
        ],
        'Reverse': [
          {'name': 'region2', 'start': 21, 'end': 58},
          {'name': 'region1', 'start': 181, 'end': 232},
        ],
      },
    },
    'awt': {
      'KL': {
        'Forward': [
          {'name': 'exon1', 'start': 21, 'end': 66},
          {'name': 'intron1', 'start': 67, 'end': 111},
          {'name': 'intron2', 'start': 139, 'end': 195},
          {'name': 'splice', 'start': 196, 'end': 201},
          {'name': 'exon2', 'start': 202, 'end': 233},
        ],
        'Forward': [
          {'name': 'exon2', 'start': 21, 'end': 52},
          {'name': 'splice', 'start': 53, 'end': 58},
          {'name': 'intron2', 'start': 59, 'end': 115},
          {'name': 'intron1', 'start': 143, 'end': 187},
          {'name': 'exon1', 'start': 188, 'end': 233},
        ],
      }
    },
    'd5': {
      'KL': {
        'Forward': [
          {'name': 'exon1', 'start': 21, 'end': 66},
          {'name': 'intron1', 'start': 67, 'end': 111},
          {'name': 'intron2', 'start': 139, 'end': 195},
          {'name': 'exon2', 'start': 195, 'end': 227},
        ],
        'Forward': [
          {'name': 'exon2', 'start': 21, 'end': 52},
          {'name': 'intron2', 'start': 53, 'end': 109},
          {'name': 'intron1', 'start': 137, 'end': 181},
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
        {'name': 'exon1', 'start': 1, 'end': 252},
      ],
      'Reverse': [
        {'name': 'exon1', 'start': 1, 'end': 252},
      ],
    },
    'db': {
      'Forward': [
        {'name': 'exon1', 'start': 1, 'end': 252},
      ],
      'Reverse': [
        {'name': 'exon1', 'start': 1, 'end': 252},
      ],
    },
    'awt': {
      'Forward': [
        {'name': 'exon1', 'start': 1, 'end': 66},
        {'name': 'intron', 'start': 67, 'end': 201},
        {'name': 'exon2', 'start': 202, 'end': 253},
      ],
      'Reverse': [
        {'name': 'exon2', 'start': 1, 'end': 52},
        {'name': 'intron', 'start': 53, 'end': 187},
        {'name': 'exon1', 'start': 188, 'end': 253},
      ],
    },
    'd5': {
      'Forward': [
        {'name': 'exon1', 'start': 1, 'end': 66},
        {'name': 'intron', 'start': 67, 'end': 195},
        {'name': 'exon2', 'start': 196, 'end': 247},
      ],
      'Reverse': [
        {'name': 'exon2', 'start': 1, 'end': 52},
        {'name': 'intron', 'start': 53, 'end': 181},
        {'name': 'exon1', 'start': 182, 'end': 247},
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
    'exon1': '#ff0000',
    'exon2': '#ff0000',
    'intron': '#548235',
  },
  'BAR_TEXT': {
    'exon1': 'Exon1',
    'exon2': 'Exon2',
    'intron': 'Intron',
  },
  'AREAS_COLOR': {
    'exon1': '#ff0000',
    'exon2': '#ff0000',
    'intron': '#548235',
    'region1': '#ff0000',
    'region2': '#548235',
  },
  'AREAS_TEXT': {
    'exon1': 'Exon1',
    'exon2': 'Exon2',
    'intron': 'Intron',
    'region1': 'Region1',
    'region2': 'Region2',
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

