import json
import pandas as pd
import os

MMEJ_LIST_CSV = os.path.join('csv', 'plasmid', 'mmej_color.csv')
SCHEMES_DATA_JS = os.path.join('html', 'plasmid', 'data.js')

JS_VARS = {
  'REF_SEQ': {
    'wt': {
      'Forward': 'GCCTCTTTAAAAGCTTGACCGAGAGCAATCCCGCAGTCTTCAGTGGTGTGATGGTCGTCTATGTGTAAGTCACCAATGCACTCAACGATTAGCGACCAGCCGGAATGCTTGGGTATGTTAATATGGACTAAAGGAGGCTTTTCTGCAGGTCGACTCTAGAACCACTCTACAAAACCAAAACCAGGGTTTATAAAATTATACTGTTGCGGAAAGCTGAAACTAAAAGAAAAACCCGACTATGCTATTTTAATCATTGAAAACGAATTTATTTAGATCCCCGTACAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTCCAACAGCCAGAGCATGTATCATATGGTCCAGAAACCCTATACCTGTGTGGACGTTAATCACTTGCGATTGTGTGGCCTGTTCTGCTACTG',
      'Reverse': 'CAGTAGCAGAACAGGCCACACAATCGCAAGTGATTAACGTCCACACAGGTATAGGGTTTCTGGACCATATGATACATGCTCTGGCTGTTGGAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTGTACGGGGATCTAAATAAATTCGTTTTCAATGATTAAAATAGCATAGTCGGGTTTTTCTTTTAGTTTCAGCTTTCCGCAACAGTATAATTTTATAAACCCTGGTTTTGGTTTTGTAGAGTGGTTCTAGAGTCGACCTGCAGAAAAGCCTCCTTTAGTCCATATTAACATACCCAAGCATTCCGGCTGGTCGCTAATCGTTGAGTGCATTGGTGACTTACACATAGACGACCATCACACCACTGAAGACTGCGGGATTGCTCTCGGTCAAGCTTTTAAAGAGGC',
    },
    'db': {
      'Forward': 'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
      'Reverse': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGTGGATAACTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
    },
  },
  'AREAS': {
    'wt': {
      '2dsb': {
        'Forward': [
          {'name': 'exon1', 'start': 21, 'end': 112},
          {'name': 'intron', 'start': 113, 'end': 287},
          {'name': 'branch', 'start': 288, 'end': 320},
          {'name': 'intron', 'start': 321, 'end': 341},
          {'name': 'exon2', 'start': 342, 'end': 405},
        ],
        'Reverse': [
          {'name': 'exon1', 'start': 21, 'end': 84},
          {'name': 'intron', 'start': 85, 'end': 105},
          {'name': 'branch', 'start': 106, 'end': 138},
          {'name': 'intron', 'start': 139, 'end': 313},
          {'name': 'exon2', 'start': 314, 'end': 405},
        ],
      }
    },
    'db': {
      '2dsb': {
        'Forward': [
          {'name': 'exon1', 'start': 21, 'end': 112},
          {'name': 'intron', 'start': 113, 'end': 308},
          {'name': 'exon2', 'start': 342, 'end': 405},
        ],
        'Reverse': [
          {'name': 'exon2', 'start': 21, 'end': 84},
          {'name': 'intron', 'start': 85, 'end': 280},
          {'name': 'exon1', 'start': 281, 'end': 372},
        ],
      }
    }
  },
  # FIXME: MAKE SURE THE NESTING IS CORRECT WITH THE 2DSB or HG39/42 KEYS
  # COMPARE WITH PLASMID DATA, THEN ROLL!!
  'PRIMER': {
    'wt': {
      'Forward': [
        {'name': 'p1', 'start': 1, 'end': 20},
        {'name': 'p2', 'start': 406, 'end': 425},
      ],
      'Reverse': [
        {'name': 'p2', 'start': 1, 'end': 20},
        {'name': 'p1', 'start': 406, 'end': 425},
      ],
    },
    'db': {
      'Forward': [
        {'name': 'p1', 'start': 1, 'end': 20},
        {'name': 'p2', 'start': 373, 'end': 392},
      ],
      'Reverse': [
        {'name': 'p2', 'start': 1, 'end': 20},
        {'name': 'p1', 'start': 373, 'end': 392},
      ],
    },
  },
  'BARS': {
    'wt': {
      'Forward': [
        {'name': 'HIS3', 'start': 21, 'end': 112},
        {'name': 'Intron', 'start': 113, 'end': 341},
        {'name': 'HIS3', 'start': 342, 'end': 405},
      ],
      'Reverse': [
        {'name': 'HIS3', 'start': 21, 'end': 84},
        {'name': 'Intron', 'start': 85, 'end': 313},
        {'name': 'HIS3', 'start': 314, 'end': 405},
      ],
    },
    'db': {
      'Forward': [
        {'name': 'HIS3', 'start': 21, 'end': 112},
        {'name': 'Intron', 'start': 113, 'end': 308},
        {'name': 'HIS3', 'start': 342, 'end': 405},
      ],
      'Reverse': [
        {'name': 'HIS3', 'start': 21, 'end': 84},
        {'name': 'Intron', 'start': 85, 'end': 280},
        {'name': 'HIS3', 'start': 281, 'end': 372},
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

