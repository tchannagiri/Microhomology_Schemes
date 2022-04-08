import json
import time
import os
from pathlib import Path
import pandas as pd
import numpy as np
import re
import seaborn as sns
import random
import itertools

MICROHOMOLOGY_TYPES = ['normal', 'pseudo', 'random']
DATA_DIR = {
  'normal': 'microhomologies_normal',
  'random': 'microhomologies_random',
  'pseudo': 'microhomologies_pseudo',
}

NUM_REPEATS = 4

WINDOW_RADIUS = 10

NUM_PSEUDO = 1000

ALL_FILES = {
  'normal': {},
  'pseudo': {},
  'random': {},
}

def init_all_files():
  for micro_type in ALL_FILES:
    ALL_FILES[micro_type]['all_microhomologies_csv'] = (
      os.path.join(DATA_DIR[micro_type], f'all_microhomologies_{micro_type}.csv')
    )
    ALL_FILES[micro_type]['all_microhomologies_txt'] = (
      os.path.join(DATA_DIR[micro_type], f'all_microhomologies_{micro_type}.txt')
    )
    ALL_FILES[micro_type]['microhomology_vars_javascript'] = (
      os.path.join('microhomology_svg_schemes', f'microhomology_vars_{micro_type}.js')
    )
    ALL_FILES[micro_type]['all_microhomologies_counts_csv'] = (
      os.path.join(DATA_DIR[micro_type], f'all_microhomologies_counts_{micro_type}.csv')
    )

LENGTH_LIST = [3, 4, 5, 6]
HGUIDES = ['HG39', 'HG42']
STRANDS = ['R1', 'R2']

MICROHOMOLOGY_GROUPS_FOR_TREATMENT = {
  'normal': {
    'spliced': ['G1', 'G2', 'G3', 'G4', 'G5'],
    'noCMV': ['G1', 'G2', 'G3', 'G4', 'G5'],
    'nonspliced': ['G1', 'G3', 'G5'],
  },
  'pseudo': {
    'spliced': ['G1X'],
    'noCMV': ['G1X'],
    'nonspliced': ['G1X'],
  },
  'random': {
    'spliced': ['GX'],
    'noCMV': ['GX'],
    'nonspliced': ['GX'],
  },
}

MICROHOMOLOGY_GROUPS_FOR_TREATMENT_PAIR = {
  'normal': {
    ('spliced', 'noCMV'): ['G1', 'G2', 'G3', 'G4', 'G5'],
    ('spliced', 'nonspliced'): ['G1', 'G3', 'G5'],
  },
  'pseudo': {
    ('spliced', 'noCMV'): ['G1X'],
    ('spliced', 'nonspliced'): ['G1X'],
  },
  'random': {
    ('spliced', 'noCMV'): ['GX'],
    ('spliced', 'nonspliced'): ['GX'],
  },
}

CELL_LINES = ['KO', 'WT']
TREATMENTS = ['spliced', 'nonspliced', 'noCMV']

TREATMENTS_FOR_MICROHOMOLOGY_GROUP = {
  'G1': ['spliced', 'noCMV', 'nonspliced'],
  'G2': ['spliced', 'noCMV'],
  'G3': ['spliced', 'noCMV', 'nonspliced'],
  'G4': ['spliced', 'noCMV'],
  'G5': ['spliced', 'noCMV', 'nonspliced'],
  'G1X': ['spliced', 'noCMV', 'nonspliced'],
  'GX': ['spliced', 'noCMV', 'nonspliced'],
}

MICROHOMOLOGY_FIELDS = [
  'microhomology_id',
  'hguide',
  'strand',
  'treatment',
  'microhomology_group',
  'match_length',
  'microhomology_match',
  'pos_left',
  'pos_right',
  'microhomology_repaired',
  'window_around_repair',
  'window_radius',
  'ref_seq',
  'visual_comparison',
  'scheme',
]

REF_SEQ = {
  'R1': 'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
  'R2': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
}

REF_SEQ_LENGTH = len(REF_SEQ['R1'])

if REF_SEQ_LENGTH != 229:
  raise ValueError(f'Wrong reference sequence length: {REF_SEQ_LENGTH}')

DISJOINT_MICROHOMOLOGY_GROUPS = {
  'normal': {
    'spliced': ['G1', 'G2'],
    'noCMV': ['G1', 'G2'],
    'nonspliced': ['G1', 'G3'],
  },
}

AREAS = {
  'HG39': {
    'R1': {
      'exon1': [(21, 67)],
      'exon2': [(184, 209)],
      'intron': [(73, 183)],
      'branch': [(119, 173)],
      'intronNoBranch': [(73, 118), (174, 183)],
      'branchFlank1': [(21, 118)],
      'branchFlank2': [(174, 209)],
      'none': [(0,0)],
    },
    'R2': None,
  },
  'HG42': {
    'R2': {
      'exon2': [(21, 46)],
      'exon1': [(158, 209)],
      'intron': [(47, 157)],
      'branch': [(57, 111)],
      'intronNoBranch': [(47, 56), (112, 157)],
      'branchFlank2': [(21, 56)],
      'branchFlank1': [(112, 209)],
      'none': [(0,0)],
    },
    'R1': None,
  }
}

MICROHOMOLOGY_AREA_SYMBOLS = {
  'G1': ('E', 'E'),
  'G2': ('E', 'I'),
  'G3': ('E', 'I'),
  'G4': ('E', 'B'),
  'G5': ('E', 'E'),
  'G1X': ('X', 'X'),
}

MICROHOMOLOGY_GROUP_PREFIXES = {
  'normal': [
    'G1',
    'G2',
    'G3',
    'G4',
    'G5',
  ],
  'pseudo': [
    'G1X',
  ],
  'random': [
    'GX',
  ],
}

MICROHOMOLOGY_GROUP_NUMBER = {
  'G1': 1,
  'G2': 2,
  'G3': 3,
  'G4': 4,
  'G5': 5,
  'GX': 6,
}

MICROHOMOLOGY_GROUP_PALETTE = sns.color_palette('deep').as_hex()

MICROHOMOLOGY_GROUP_COLOR = dict(zip(
  (
      MICROHOMOLOGY_GROUP_PREFIXES['normal'] +
      MICROHOMOLOGY_GROUP_PREFIXES['pseudo'] +
      MICROHOMOLOGY_GROUP_PREFIXES['random']
  ),
  MICROHOMOLOGY_GROUP_PALETTE,
))

TREATMENT_PALETTE = sns.color_palette('muted').as_hex()

TREATMENT_COLOR = dict(zip(TREATMENTS, TREATMENT_PALETTE))

MICROHOMOLOGY_GROUP_DESCRIPTION = {
  'G1': 'Exon / Exon',
  'G2': 'Exon / Intron (All)',
  'G3': 'Exon / Intron (No Branch-Site)',
  'G4': 'Exon / Branch-Site',
  'G5': 'Exon / Exon (No 6bp Microhomology)',
  'GX': 'Random Sequences',
  'G1X': 'Branch-Flank-1 / Branch-Frank-2',
}

MICROHOMOLOGY_GROUP_NAMES = {
  'HG39': {
    'R1': { 'G1': None, 'G2': None, 'G3': None, 'G4': None, 'G5': None, 'G1X': None},
    'R2': { 'G1': None, 'G2': None, 'G3': None, 'G4': None, 'G5': None, 'G1X': None},
  },
  'HG42': {
    'R1': { 'G1': None, 'G2': None, 'G3': None, 'G4': None, 'G5': None, 'G1X': None},
    'R2': { 'G1': None, 'G2': None, 'G3': None, 'G4': None, 'G5': None, 'G1X': None},
  },
}

MICROHOMOLOGY_GROUP_AREA_PAIRS = {
  'HG39': {
    'R1': {
      'G1': {'left': None, 'right': None},
      'G2': {'left': None, 'right': None},
      'G3': {'left': None, 'right': None},
      'G4': {'left': None, 'right': None},
      'G5': {'left': None, 'right': None},
      'G1X': {'left': None, 'right': None},
      'GX': {'left': None, 'right': None},
    },
    'R2': {
      'G1': {'left': None, 'right': None},
      'G2': {'left': None, 'right': None},
      'G3': {'left': None, 'right': None},
      'G4': {'left': None, 'right': None},
      'G5': {'left': None, 'right': None},
      'G1X': {'left': None, 'right': None},
      'GX': {'left': None, 'right': None},
    },
  },
  'HG42': {
    'R1': {
      'G1': {'left': None, 'right': None},
      'G2': {'left': None, 'right': None},
      'G3': {'left': None, 'right': None},
      'G4': {'left': None, 'right': None},
      'G5': {'left': None, 'right': None},
      'G1X': {'left': None, 'right': None},
      'GX': {'left': None, 'right': None},
    },
    'R2': {
      'G1': {'left': None, 'right': None},
      'G2': {'left': None, 'right': None},
      'G3': {'left': None, 'right': None},
      'G4': {'left': None, 'right': None},
      'G5': {'left': None, 'right': None},
      'G1X': {'left': None, 'right': None},
      'GX': {'left': None, 'right': None},
    },
  },
}

MICROHOMOLOGY_GROUP_POST_FILTER = {
  'G1': None,
  'G2': None,
  'G3': None,
  'G4': None,
  'G5': lambda micro: len(micro['microhomology_match']) < 6,
  'G1X': None,
  'GX': None,
}

PSEUDO_GROUP = {
  'G1': 'normal',
  'G2': 'normal',
  'G3': 'normal',
  'G4': 'normal',
  'G5': 'normal',
  'GX': 'pseudo',
}

TREATMENT_DELETE_BRANCH = {
  'spliced': False,
  'noCMV': False,
  'nonspliced': True,
}

HGUIDE_ISOLATED_EXON = {
  'HG39': 'exon1',
  'HG42': 'exon2',
}

STRAND_EXON_ORDER = {
  'R1': ('exon1', 'exon2'),
  'R2': ('exon2', 'exon1'),
}

RANDOM_RADIUS = 13

def make_arg_str(*args_1, **args_2):
  arg_str_1 = ''.join([f'_[{x}]' for x in args_1 if x is not None])
  arg_str_2 = ''.join([f'_[{key}={value}]' for key, value in args_2.items() if value is not None])
  return arg_str_1 + arg_str_2

def reverse_complement(dna):
  complement_letter = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
  }
  dna = list(dna)[::-1]
  dna = [complement_letter[x] for x in dna]
  return ''.join(dna)

def init_R2_ref_seq():
  REF_SEQ['R2'] = reverse_complement(REF_SEQ['R1'])

def init_reverse_strand_areas():
  for hguide in AREAS:
    old_strand = list(AREAS[hguide])[0]
    new_strand = 'R2' if old_strand == 'R1' else 'R1'
    AREAS[hguide][new_strand] = {}
    for area_name in AREAS[hguide][old_strand]:
      new_area_list = []
      for old_area in AREAS[hguide][old_strand][area_name]:
        new_area = (
          REF_SEQ_LENGTH - old_area[1] + 1,
          REF_SEQ_LENGTH - old_area[0] + 1,
        )
        new_area_list.append(new_area)
      AREAS[hguide][new_strand][area_name] = new_area_list

def init_microhomology_group_area_pairs():
  for micro_type in MICROHOMOLOGY_TYPES:
    for hguide in HGUIDES:
      for strand in STRANDS:
        for group_prefix in MICROHOMOLOGY_GROUP_PREFIXES[micro_type]:
          if group_prefix in ['G1', 'G5']:
            left = STRAND_EXON_ORDER[strand][0]
            right = STRAND_EXON_ORDER[strand][1]
          elif group_prefix == 'G2':
            left = HGUIDE_ISOLATED_EXON[hguide]
            right = 'intron'
          elif group_prefix == 'G3':
            left = HGUIDE_ISOLATED_EXON[hguide]
            right = 'intronNoBranch'
          elif group_prefix == 'G4':
            left = HGUIDE_ISOLATED_EXON[hguide]
            right = 'branch'
          elif group_prefix == 'G1X':
            left = 'branchFlank1'
            right = 'branchFlank2'
          elif group_prefix == 'GX':
            left = 'none'
            right = 'none'
          if AREAS[hguide][strand][left][0][0] > AREAS[hguide][strand][right][0][0]:
            left, right = right, left
          MICROHOMOLOGY_GROUP_AREA_PAIRS[hguide][strand][group_prefix]['left'] = left
          MICROHOMOLOGY_GROUP_AREA_PAIRS[hguide][strand][group_prefix]['right'] = right
      
def init_microhomology_group_names():
  for micro_type in MICROHOMOLOGY_TYPES:
    for hguide in HGUIDES:
      for strand in STRANDS:
        for group_prefix in MICROHOMOLOGY_GROUP_PREFIXES[micro_type]:
          if group_prefix == 'GX':
            MICROHOMOLOGY_GROUP_NAMES[hguide][strand][group_prefix] = '_'.join([
              hguide,
              strand,
              group_prefix,
              'random',
            ])
          else:
            MICROHOMOLOGY_GROUP_NAMES[hguide][strand][group_prefix] = '_'.join([
              hguide,
              strand,
              group_prefix,
              MICROHOMOLOGY_GROUP_AREA_PAIRS[hguide][strand][group_prefix]['left'],
              MICROHOMOLOGY_GROUP_AREA_PAIRS[hguide][strand][group_prefix]['right'],
            ])

def log(*strs):
  print(': '.join([time.strftime("%H:%M:%S", time.localtime()), *strs]))

def create_parent_dir(file_name):
   Path(os.path.dirname(file_name)).mkdir(parents=True, exist_ok=True)

def write_csv(data: pd.DataFrame, file_name):
  create_parent_dir(file_name)
  data.to_csv(file_name, index=False)

def read_csv(file_name):
  return pd.read_csv(file_name)

def extract_window(ref_seq, pos_left, pos_right, match_length, delete_branch_site, window_radius):
  pos_left_0 = pos_left - 1
  pos_right_0 = pos_right - 1

  if delete_branch_site is not None:
    need_delete_branch_site = pos_right <= delete_branch_site[0]
  else:
    need_delete_branch_site = False

  if need_delete_branch_site:
    # The branch site is not deleted by default but we must force it to be
    branch_site_start_0 = delete_branch_site[0] - 1
    branch_site_end_0 = delete_branch_site[1] - 1
    repaired_seq = (
      ref_seq[: pos_left_0] +
      ref_seq[pos_right_0 : branch_site_start_0] +
      ref_seq[branch_site_end_0 + 1 :]
    )
    visual_comparison = (
      ('+' * pos_left_0) +
      ('M' * match_length) +
      ('-' * (pos_right_0 - pos_left_0 - match_length)) +
      ('M' * match_length) +
      ('+' * (branch_site_start_0 - pos_right_0 - match_length)) +
      ('-' * (branch_site_end_0 - branch_site_start_0 + 1)) +
      ('+' * (len(ref_seq) - branch_site_end_0 - 1))
    )
  else:
    repaired_seq = ref_seq[:pos_left_0] + ref_seq[pos_right_0:]
    visual_comparison = (
      ('+' * pos_left_0) +
      ('M' * match_length) +
      ('-' * (pos_right_0 - pos_left_0 - match_length)) +
      ('M' * match_length) +
      ('+' * (len(ref_seq) - pos_right_0 - match_length))
    )

  window_start = max(pos_left_0 - window_radius, 0)
  window_end = min(pos_left_0 + window_radius + match_length, len(repaired_seq))
  window = repaired_seq[window_start : window_end]

  micro_data = {field: None for field in MICROHOMOLOGY_FIELDS}
  micro_data.update({
    'match_length': match_length,
    'microhomology_match': ref_seq[pos_left_0 : pos_left_0 + match_length],
    'pos_left': pos_left,
    'pos_right': pos_right,
    'microhomology_repaired': repaired_seq,
    'window_around_repair': window,
    'window_radius': window_radius,
    'ref_seq': ref_seq,
    'visual_comparison': visual_comparison,
  })
  return micro_data

# Finds the matching substring in area_1 and area_2 of a given length
# and returns the two start positions, the microhomology match string, the
# repaired sequence and a visualization of how the reference sequence is repaired
# (ie. which portions are deleted from the reference)
# areas_1 and areas_2 are given in 1-based coordinates
# returned positions are also in 1-based coordinates
def make_microhomologies(
  ref_seq,
  area_left,
  area_right,
  length_list,
  delete_branch_site,
  window_radius,
):
  microhomologies = []
  start_left, end_left = area_left
  start_right, end_right = area_right
  for match_length in length_list:
    for pos_left in range(start_left, end_left + 1):
      for pos_right in range(start_right, end_right + 1):
        if (pos_left + match_length - 1 > end_left) or \
          (pos_right + match_length - 1 > end_right):
          continue
        pos_left_0 = pos_left - 1
        pos_right_0 = pos_right - 1
        match_1 = ref_seq[pos_left_0 : pos_left_0 + match_length]
        match_2 = ref_seq[pos_right_0 : pos_right_0 + match_length]
        if match_1 == match_2:
          extendable_left = (
            (pos_left - 1 >= start_left) and
            (pos_right - 1 >= start_right) and
            (ref_seq[pos_left_0 - 1] == ref_seq[pos_right_0 - 1])
          )
          extendable_right = (
            (pos_left + match_length <= end_left) and
            (pos_right + match_length <= end_right) and
            (ref_seq[pos_left_0 + match_length] == ref_seq[pos_right_0 + match_length])
          )
          if (not extendable_left) and (not extendable_right):
            microhomologies.append(extract_window(
              ref_seq,
              pos_left,
              pos_right,
              match_length,
              delete_branch_site,
              window_radius,
            ))
  return microhomologies

# Extract all windows in the left and right
# def make_pseudo_microhomologies_all_flank(
#   ref_seq,
#   area_left,
#   area_right,
#   delete_branch_site,
#   window_radius,
# ):
#   microhomologies = []
#   start_left, end_left = area_left
#   start_right, end_right = area_right
#   for pos_left in range(start_left, end_left + 1):
#     for pos_right in range(start_right, end_right + 1):
#       micro_data = extract_window(
#         ref_seq,
#         pos_left,
#         pos_right,
#         0,
#         delete_branch_site,
#         window_radius,
#       )
#       micro_data['microhomology_match'] = 'X'
#       microhomologies.append(micro_data)
#       break
#     break
#   return microhomologies

def sample_pos_pairs(n, area_left, area_right, seed=0):
  random.seed(seed)

  num_possible = (area_left[1] - area_left[0] + 1) * (area_right[1] - area_right[0] + 1)
  if n > num_possible:
    raise ValueError(f'Cannot sample {n} values out of a population size {num_possible}.')
  elif n > num_possible // 2:
    all_possible = list(itertools.product(
      range(area_left[0], area_left[1] + 1),
      range(area_right[0], area_right[1] + 1),
    ))
    return random.sample(all_possible, n)
  else:
    samples = set()
    while len(samples) < n:
      left_pos = random.randrange(area_left[0], area_left[1] + 1)
      right_pos = random.randrange(area_right[0], area_right[1] + 1)
      samples.add((left_pos, right_pos))
    return list(samples)

# Extract 100 random windows by simulating the repair of right and left
# avoid the normal microhomologies
def make_pseudo_microhomologies(
  normal_microhomologies,
  ref_seq,
  area_left,
  area_right,
  delete_branch_site,
  window_radius,
  size = None,
):
  microhomologies = []
  if size is None:
    pos_sample = list(itertools.product(
      range(area_left[0], area_left[1] + 1),
      range(area_right[0], area_right[1] + 1),
    ))
  else:
    pos_sample = sample_pos_pairs(size, area_left, area_right, seed=area_left[0])
  # pos_normal = {
  #   (x['pos_left'] + d, x['pos_right'] + d)
  #   for x in normal_microhomologies
  #   for d in range(x['match_length'] + 1)
  # }
  # pos_sample = set(pos_sample) - pos_normal
  for pos_left, pos_right in pos_sample:
    micro_data = extract_window(
      ref_seq,
      pos_left,
      pos_right,
      0,
      delete_branch_site,
      window_radius,
    )
    micro_data['microhomology_match'] = 'X'
    microhomologies.append(micro_data)
  microhomologies = [
    x
    for x in microhomologies
    if not any(
      (x['window_around_repair'] in y['window_around_repair']) or
      (y['window_around_repair'] in x['window_around_repair'])
      for y in normal_microhomologies
    )
  ]
  return microhomologies

# Just get random repair windows
def make_random_microhomologies(window_radius):
  microhomologies = []
  DNA_letters = 'ACGT'
  np.random.seed(0)
  for _ in range(1000):
    random_ints = np.random.randint(0, 4, size=(window_radius * 2))
    repair_window = ''.join([
      DNA_letters[x] for x in random_ints
    ])
    micro_data = {field: None for field in MICROHOMOLOGY_FIELDS}
    micro_data.update({
      'window_around_repair': repair_window,
      'window_radius': window_radius,
    })
    microhomologies.append(micro_data)
  return microhomologies

def pretty_print_microhomologies(microhomologies, file_name_out):
  field_list = [
    x for x in MICROHOMOLOGY_FIELDS
    if x not in ['ref_seq', 'visual_comparison', 'scheme']
  ]
  create_parent_dir(file_name_out)
  with open(file_name_out, 'w') as file_out:
    for data in microhomologies.to_dict('records'):
      for field in field_list:
        file_out.write(f'{field}: {data[field]}\n')
      file_out.write('visual_comparison:\n')
      file_out.write(data['ref_seq'] + '\n')
      file_out.write(data['visual_comparison'] + '\n')
      file_out.write(data['scheme'] + '\n')
      file_out.write('\n')

def print_javascript_microhomologies(microhomologies, file_out_name, var_name):
  with open(file_out_name, 'w') as file_out:
    file_out.write(
      f'var {var_name} = {json.dumps(microhomologies.to_dict("records"), indent=2)};\n\n'
    )

def print_javascript_extra_vars(file_out_name):
  with open(file_out_name, 'a') as file_out:
    file_out.write(f'var AREAS = {json.dumps(AREAS, indent=2)};\n\n')
    file_out.write(f'var MICROHOMOLOGY_GROUP_NAMES = {json.dumps(MICROHOMOLOGY_GROUP_NAMES, indent=2)};\n\n')
    file_out.write(f'var REF_SEQ = {json.dumps(REF_SEQ, indent=2)};\n\n')

def get_microhomology_file(pseudo, treatment, microhomology_group_name, ext):
  file_main = f'microhomologies_{treatment}_{microhomology_group_name}'
  return os.path.join(DATA_DIR[pseudo], ext, file_main + os.path.extsep + ext)

def load_all_microhomologies(pseudo):
  microhomologies = pd.read_csv(ALL_FILES[pseudo]['all_microhomologies_csv'])
  return microhomologies

def get_microhomology_id(
  treatment,
  microhomology_group_name,
  pos_1,
  pos_2,
  microhomology_match,
):
  return '_'.join([
    f'{treatment}',
    microhomology_group_name,
    f'P{pos_1}_P{pos_2}',
    microhomology_match,
  ])

def get_random_microhomology_id(
  treatment,
  microhomology_group_name,
  microhomology_idx,
):
  return '_'.join([
    f'{treatment}',
    microhomology_group_name,
    f'N{microhomology_idx}',
  ])

def remove_treatment_from_microhomology_id(microhomology_id: str):
  for treatment in TREATMENTS:
    microhomology_id = microhomology_id.replace(treatment + '_', '')
  return microhomology_id

def get_pos_from_microhomology_id(microhomology_id):
  pos = filter(lambda x: x[0] == 'P', microhomology_id.split('_'))
  pos = map(lambda x: int(x[1:]))
  return list(pos)

def make_text_scheme(length, left_symbol, left_areas, right_symbol, right_areas):
  left_areas = sorted(list(set(left_areas)), key=lambda x: x[0])
  right_areas = sorted(list(set(right_areas)), key=lambda x: x[0])
  all_areas = (
    [{'area': area, 'name': left_symbol} for area in left_areas] +
    [{'area': area, 'name': right_symbol} for area in right_areas]
  )
  all_areas_iter = iter(all_areas)
  curr_area = next(all_areas_iter, None)
  scheme_str = ''
  for i in range(1, length + 1):
    if (curr_area is not None) and (i > curr_area['area'][1]):
      curr_area = next(all_areas_iter, None)
    if (curr_area is None) or (i < curr_area['area'][0]):
      scheme_str += ' '
    elif i <= curr_area['area'][1]:
      scheme_str += curr_area['name']
    else:
      raise Exception('Impossible value for i: ' + str(i) + ', curr_area: ' + str(curr_area))
  return scheme_str

def reverse_complement(dna):
  complement_letter = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
  }
  dna = list(dna)[::-1]
  dna = [complement_letter[x] for x in dna]
  return ''.join(dna)

def make_all_microhomologies():
  all_microhomologies = {}
  for micro_type in MICROHOMOLOGY_TYPES:
    all_microhomologies[micro_type] = []
    for hguide in HGUIDES:
      for strand in STRANDS:
        ref_seq = REF_SEQ[strand]
        for group_prefix in MICROHOMOLOGY_GROUP_PREFIXES[micro_type]:
          group_name = MICROHOMOLOGY_GROUP_NAMES[hguide][strand][group_prefix]
          left_areas = AREAS[hguide][strand][
            MICROHOMOLOGY_GROUP_AREA_PAIRS[hguide][strand][group_prefix]['left']
          ]
          right_areas = AREAS[hguide][strand][
            MICROHOMOLOGY_GROUP_AREA_PAIRS[hguide][strand][group_prefix]['right']
          ]
          if micro_type == 'random':
            group_scheme = None
          else:
            group_scheme = make_text_scheme(
              len(ref_seq),
              MICROHOMOLOGY_AREA_SYMBOLS[group_prefix][0],
              left_areas,
              MICROHOMOLOGY_AREA_SYMBOLS[group_prefix][1],
              right_areas,
            )

          for treatment in TREATMENTS_FOR_MICROHOMOLOGY_GROUP[group_prefix]:
            if TREATMENT_DELETE_BRANCH[treatment]:
              delete_branch_site = AREAS[hguide][strand]['branch'][0]
            else:
              delete_branch_site = None
            for area_left in left_areas:
              for area_right in right_areas:
                if micro_type == 'random':
                  microhomologies = make_random_microhomologies(RANDOM_RADIUS)
                elif micro_type == 'pseudo':
                  normal_microhomologies = [
                    x for x in all_microhomologies['normal']
                    if x['hguide'] == hguide and
                    x['strand'] == strand and
                    x['treatment'] == treatment
                  ]
                  microhomologies = make_pseudo_microhomologies(
                    normal_microhomologies,
                    ref_seq,
                    area_left,
                    area_right,
                    delete_branch_site,
                    WINDOW_RADIUS,
                    size = NUM_PSEUDO,
                  )
                elif micro_type == 'normal':
                  microhomologies = make_microhomologies(
                    ref_seq,
                    area_left,
                    area_right,
                    LENGTH_LIST,
                    delete_branch_site,
                    WINDOW_RADIUS,
                  )
                else:
                  raise ValueError('Impossible: ' + str(micro_type))

                if MICROHOMOLOGY_GROUP_POST_FILTER[group_prefix]:
                  microhomologies = list(filter(
                    MICROHOMOLOGY_GROUP_POST_FILTER[group_prefix],
                    microhomologies
                  ))

                for micro_idx, micro in enumerate(microhomologies):
                  if micro_type == 'random':
                    micro['microhomology_id'] = get_random_microhomology_id(
                      treatment,
                      group_name,
                      micro_idx,
                    )
                  elif micro_type in ['normal', 'pseudo']:
                    micro['microhomology_id'] = get_microhomology_id(
                      treatment,
                      group_name,
                      micro['pos_left'],
                      micro['pos_right'],
                      micro['microhomology_match'],
                    )
                  else:
                    raise ValueError('Impossible: ' + str(micro_type))
                  micro['hguide'] = hguide
                  micro['strand'] = strand
                  micro['treatment'] = treatment
                  micro['microhomology_group'] = group_prefix
                  if micro_type in ['normal', 'pseudo']:
                    micro['scheme'] = group_scheme
                
                microhomologies = pd.DataFrame.from_records(microhomologies)
                if micro_type in ['normal', 'pseudo']:
                  microhomologies = microhomologies.sort_values(['pos_left', 'pos_right'])
                write_csv(microhomologies, get_microhomology_file(micro_type, treatment, group_name, 'csv'))
                if micro_type in ['normal', 'pseudo']:
                  pretty_print_microhomologies(
                    microhomologies,
                    get_microhomology_file(micro_type, treatment, group_name, 'txt'),
                  )
                all_microhomologies[micro_type] += microhomologies.to_dict('records')

  for micro_type in MICROHOMOLOGY_TYPES:
    all_microhomologies[micro_type] = pd.DataFrame.from_records(all_microhomologies[micro_type])
    write_csv(
      all_microhomologies[micro_type],
      ALL_FILES[micro_type]['all_microhomologies_csv'],
    )
    if micro_type in ['normal', 'pseudo']:
      pretty_print_microhomologies(
        all_microhomologies[micro_type],
        ALL_FILES[micro_type]['all_microhomologies_txt']
      )

      print_javascript_microhomologies(
        all_microhomologies[micro_type],
        ALL_FILES[micro_type]['microhomology_vars_javascript'],
        'MICROHOMOLOGIES',
      )
      print_javascript_extra_vars(ALL_FILES[micro_type]['microhomology_vars_javascript'])

      summary = all_microhomologies[micro_type][[
        'hguide',
        'strand',
        'treatment',
        'microhomology_group',
        'match_length',
      ]].copy()
      summary['count'] = 1
      summary = summary.groupby(
        ['hguide', 'strand', 'treatment', 'microhomology_group', 'match_length']
      ).count().reset_index()
      write_csv(summary, ALL_FILES[micro_type]['all_microhomologies_counts_csv'])

def init_all():
  init_all_files()
  init_R2_ref_seq()
  init_reverse_strand_areas()
  init_microhomology_group_area_pairs()
  init_microhomology_group_names()

init_all()

if __name__ == '__main__':
  make_all_microhomologies()

# Reference sequence for R1 files
# 5’- TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC
# Reference sequence for R2 files
# 5’- GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA