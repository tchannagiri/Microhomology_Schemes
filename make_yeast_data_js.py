import json
import pandas as pd
import os

MMEJ_LIST_CSV = os.path.join('csv', 'yeast', 'mmej_color.csv')
SCHEMES_DATA_JS = os.path.join('html', 'yeast', 'data.js')

JS_VARS = {
  'BOLDED': {
    'WT': {
      'wt': {
        '2dsb': ['EE58', 'EE63', 'EE64', 'EE65', 'EE67', 'EE58R', 'EE63R', 'EE64R', 'EE65R']
      },
      'db': {
        '2dsb': ['EE58', 'EE63', 'EE64', 'EE65', 'EE67', 'EE58R', 'EE63R', 'EE64R', 'EE65R']
      },
    },
    'KO': {
      'wt': {
        '2dsb': ['EE63', 'EE63R', 'EE65R']
      },
      'db': {
        '2dsb': ['EE63', 'EE63R', 'EE65R']
      },
    }
  },
  'REF_SEQ': {
    'wt': {
      'Forward': 'GCCTCTTTAAAAGCTTGACCGAGAGCAATCCCGCAGTCTTCAGTGGTGTGATGGTCGTCTATGTGTAAGTCACCAATGCACTCAACGATTAGCGACCAGCCGGAATGCTTGGGTATGTTAATATGGACTAAAGGAGGCTTTTCTGCAGGTCGACTCTAGAACCACTCTACAAAACCAAAACCAGGGTTTATAAAATTATACTGTTGCGGAAAGCTGAAACTAAAAGAAAAACCCGACTATGCTATTTTAATCATTGAAAACGAATTTATTTAGATCCCCGTACAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTCCAACAGCCAGAGCATGTATCATATGGTCCAGAAACCCTATACCTGTGTGGACGTTAATCACTTGCGATTGTGTGGCCTGTTCTGCTACTG',
      'Reverse': 'CAGTAGCAGAACAGGCCACACAATCGCAAGTGATTAACGTCCACACAGGTATAGGGTTTCTGGACCATATGATACATGCTCTGGCTGTTGGAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTGTACGGGGATCTAAATAAATTCGTTTTCAATGATTAAAATAGCATAGTCGGGTTTTTCTTTTAGTTTCAGCTTTCCGCAACAGTATAATTTTATAAACCCTGGTTTTGGTTTTGTAGAGTGGTTCTAGAGTCGACCTGCAGAAAAGCCTCCTTTAGTCCATATTAACATACCCAAGCATTCCGGCTGGTCGCTAATCGTTGAGTGCATTGGTGACTTACACATAGACGACCATCACACCACTGAAGACTGCGGGATTGCTCTCGGTCAAGCTTTTAAAGAGGC',
    },
    'db': {
      'Forward': 'GCCTCTTTAAAAGCTTGACCGAGAGCAATCCCGCAGTCTTCAGTGGTGTGATGGTCGTCTATGTGTAAGTCACCAATGCACTCAACGATTAGCGACCAGCCGGAATGCTTGGGTATGTTAATATGGACTAAAGGAGGCTTTTCTGCAGGTCGACTCTAGAACCACTCTACAAAACCAAAACCAGGGTTTATAAAATTATACTGTTGCGGAAAGCTGAAACTAAAAGAAAAACCCGACTATGCTATTTTAATCATTGAAAACGAATTTATTTAGATCCCCGTACAGGAAAATGGTATTATTTCCAACAGCCAGAGCATGTATCATATGGTCCAGAAACCCTATACCTGTGTGGACGTTAATCACTTGCGATTGTGTGGCCTGTTCTGCTACTG',
      'Reverse': 'CAGTAGCAGAACAGGCCACACAATCGCAAGTGATTAACGTCCACACAGGTATAGGGTTTCTGGACCATATGATACATGCTCTGGCTGTTGGAAATAATACCATTTTCCTGTACGGGGATCTAAATAAATTCGTTTTCAATGATTAAAATAGCATAGTCGGGTTTTTCTTTTAGTTTCAGCTTTCCGCAACAGTATAATTTTATAAACCCTGGTTTTGGTTTTGTAGAGTGGTTCTAGAGTCGACCTGCAGAAAAGCCTCCTTTAGTCCATATTAACATACCCAAGCATTCCGGCTGGTCGCTAATCGTTGAGTGCATTGGTGACTTACACATAGACGACCATCACACCACTGAAGACTGCGGGATTGCTCTCGGTCAAGCTTTTAAAGAGGC',
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
          {'name': 'exon2', 'start': 309, 'end': 372},
        ],
        'Reverse': [
          {'name': 'exon2', 'start': 21, 'end': 84},
          {'name': 'intron', 'start': 85, 'end': 280},
          {'name': 'exon1', 'start': 281, 'end': 372},
        ],
      }
    }
  },
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
        {'name': 'his3', 'start': 1, 'end': 112},
        {'name': 'intron', 'start': 113, 'end': 341},
        {'name': 'his3', 'start': 342, 'end': 425},
      ],
      'Reverse': [
        {'name': 'his3', 'start': 1, 'end': 84},
        {'name': 'intron', 'start': 85, 'end': 313},
        {'name': 'his3', 'start': 314, 'end': 425},
      ],
    },
    'db': {
      'Forward': [
        {'name': 'his3', 'start': 1, 'end': 112},
        {'name': 'intron', 'start': 113, 'end': 308},
        {'name': 'his3', 'start': 309, 'end': 392},
      ],
      'Reverse': [
        {'name': 'his3', 'start': 1, 'end': 84},
        {'name': 'intron', 'start': 85, 'end': 280},
        {'name': 'his3', 'start': 281, 'end': 392},
      ],
    },
  },
  'CUT_POS': {
    'wt': {
      '2dsb': {
        'Forward': [105, 340], 
        'Reverse':  [85, 320],
      },
    },
    'db': {
      '2dsb': {
        'Forward': [105, 307],
        'Reverse': [85, 287],
      },
    },
  },
  'BAR_COLOR': {
    'his3': '#4472c4',
    'intron': '#00b050',
  },
  'BAR_TEXT': {
    'his3': 'HIS3',
    'intron': 'Intron',
  },
  'AREAS_COLOR': {
    'exon1': '#4472c4', 
    'exon2': '#4472c4',
    'intron': '#00b050',
    'branch': '#548235',
  },
  'AREAS_TEXT': {
    'exon1': 'Exon1',
    'exon2': 'Exon2',
    'intron': 'Intron',
    'branch': 'Intron (branch site)',
  },
  'PRIMER_TEXT': {
    'p1': 'Primer',
    'p2': 'Primer',
  },
  'LABELS': {
    '2dsb': '2 DSB',
    'Forward': 'forward strand',
    'Reverse': 'reverse strand',
    'wt': 'Antisense',
    'db': 'BranchΔ',
    'all': None,
    'exon_exon': 'Exon-Exon group',
    'exon_intron': 'Exon-Intron group',
  },
  'CELLTYPE_FLIPPED_INTRON': ['wt', 'db'],
}

microhomologies = [
  pd.read_csv(MMEJ_LIST_CSV)
]
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
