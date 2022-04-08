var PCIS = {
  'R1': [
    ['f6', 1, 20],
    ['r5', 210, 229],
  ],
  'R2': [
    ['r5', 1, 20],
    ['f6', 210, 229],
  ],
}

var BARS = {
  'R1': [
    ['dsred', 1, 72],
    ['intron', 73, 183],
    ['dsred', 184, 229],
  ],
  'R2': [
    ['dsred', 1, 46],
    ['intron', 47, 157],
    ['dsred', 158, 229],
  ],
};

var BAR_COLOR = {
  'dsred': 'red',
  'intron': 'pink',
};

var BAR_TEXT = {
  'dsred': 'DS Red',
  'intron': 'Intron',
};

var AREAS_COLOR = {
  'exon1': 'blue',
  'exon2': 'blue',
  'intron': 'red',
  'branch': 'black',
};

var AREAS_TEXT = {
  'exon1': 'Exon 1',
  'exon2': 'Exon 2',
  'intron': 'Intron',
  'branch': 'Intron (branch site)',
};

var PCIS_TEXT = {
  'r5': 'DsRed.pCis.R5',
  'f6': 'DsRed.pCis.F6',
};

var PCIS_COLOR = 'purple';

var DELETE_BRANCH = {
  'wt': false,
  'db': true,
};

var SCALE_X = 10;
var SCALE_Y = 20;

function reverse_complement_dna(s) {
  var t = '';
  for (var i = s.length - 1; i >= 0; i--) {
    switch (s[i]) {
      case 'A': t += 'T'; break;
      case 'C': t += 'G'; break;
      case 'G': t += 'C'; break;
      case 'T': t += 'A'; break;
      default: throw `Bad DNA letter: ${s[i]}`;
    }
  }
  return t;
}

function init_svgs_separate(content_selector, microhomologies, labels_width, bar_width, scheme_height, padding) {
  return (
    d3.select(content_selector)
    .selectAll('svg')
    .data(microhomologies)
    .enter()
    .append('svg')
    .attr('width', labels_width + bar_width)
    .attr('height', scheme_height + padding)
  );
}

function init_svgs_combined(content_selector, labels_width, bar_width, scheme_height, padding) {
  return (
    d3.select(content_selector)
      .append('svg')
      .attr('width', labels_width + bar_width)
      .attr('height', scheme_height + padding)
  );
}

function make_seq_svgs_separate(
  svgs,
  labels_width,
  bar_width,
  seq_y,
  seq_height,
  ref_seq,
  delete_branch,
) {
  var seq_svgs = svgs.append('svg')
    .style('overflow', 'visible')
    .attr('x', labels_width)
    .attr('y', seq_y)
    .attr('width', bar_width)
    .attr('height', seq_height);

  seq_svgs.append('text')
    .attr('x', 0)
    .attr('y', 0)
    .attr("dy", "1em")
    .attr('height', seq_height)
    .attr('textLength', bar_width)
    .style('font-family', 'monospace')
    .text(ref_seq);

  for (var pos of ['pos_left', 'pos_right']) {
    seq_svgs.append('rect')
      .attr('x', d => (d[pos]- 1) * SCALE_X)
      .attr('y', 0)
      .attr('width', d => d['microhomology_match'].length * SCALE_X)
      .attr('height', seq_height)
      .style('fill', 'transparent')
      .style('stroke', 'black')
      .style('stroke-width', 2);
  }

  seq_svgs.append('rect')
    .attr('x', d => (d['pos_left'] + d['microhomology_match'].length - 1) * SCALE_X)
    .attr('y', 0)
    .attr('width', d => (d['pos_right'] - d['pos_left'] - d['microhomology_match'].length) * SCALE_X)
    .attr('height', seq_height)
    .attr('opacity', 0.3)
    .style('fill', 'black')
    .style('stroke', 'transparent');

  if (delete_branch) {
    seq_svgs.filter(d => (d['pos_left'] > delete_branch[1]) || (d['pos_right'] < delete_branch[0]))
      .append('rect')
      .attr('x', (delete_branch[0] - 1) * SCALE_X)
      .attr('y', 0)
      .attr('width', (delete_branch[1] - delete_branch[0] + 1) * SCALE_X)
      .attr('height', seq_height)
      .attr('opacity', 0.3)
      .style('fill', 'black')
      .style('stroke', 'transparent');
  }

  return seq_svgs;
}

function make_seq_svgs_combined(
  svgs,
  labels_width,
  bar_width,
  seq_y,
  seq_height,
  ref_seq,
  microhomologies,
  delete_branch,
) {
  var label_svgs = svgs.append('svg')
    .style('overflow', 'visible')
    .attr('x', 0)
    .attr('y', seq_y)
    .attr('width', labels_width)
    .attr('height', seq_height * microhomologies.length);

  var seq_svgs = svgs.append('svg')
    .style('overflow', 'visible')
    .attr('x', labels_width)
    .attr('y', seq_y)
    .attr('width', bar_width)
    .attr('height', seq_height * microhomologies.length);

  for (var i = 0; i < microhomologies.length; i++) {
    label_svgs.append('svg')
      .attr('x', 0)
      .attr('y', i * seq_height)
      .attr('width', labels_width)
      .attr('height', seq_height)
      .style('overflow', 'visible')
      .append('text')
      .attr('width', labels_width)
      .attr('height', seq_height)
      .attr("dy", '1em')
      .style('font-family', 'monospace')
      .text(microhomologies[i]['microhomology_id']);

    seq_svgs.append('text')
      .attr('x', 0)
      .attr('y', i * seq_height)
      .attr("dy", "1em")
      .attr('height', seq_height)
      .attr('textLength', bar_width)
      .style('font-family', 'monospace')
      .text(ref_seq);
  }

  var micro_idx = 0;
  for (var micro of microhomologies) {
    for (var pos of ['pos_left', 'pos_right']) {
      seq_svgs.append('rect')
        .attr('x', (micro[pos] - 1) * SCALE_X)
        .attr('y', micro_idx * seq_height)
        .attr('width', micro['microhomology_match'].length * SCALE_X)
        .attr('height', seq_height)
        .style('fill', 'transparent')
        .style('stroke', 'black')
        .style('stroke-width', 2);
    }

    seq_svgs.append('rect')
      .attr('x', (micro['pos_left'] + micro['microhomology_match'].length - 1) * SCALE_X)
      .attr('y', micro_idx * seq_height)
      .attr('width', (micro['pos_right'] - micro['pos_left'] - micro['microhomology_match'].length) * SCALE_X)
      .attr('height', seq_height)
      .attr('opacity', 0.3)
      .style('fill', 'black')
      .style('stroke', 'transparent');

    if (
      delete_branch &&
      ((micro['pos_left'] > delete_branch[1]) || (micro['pos_right'] < delete_branch[0]))
    ) {
      seq_svgs.append('rect')
        .attr('x', (delete_branch[0] - 1) * SCALE_X)
        .attr('y', micro_idx * seq_height)
        .attr('width', (delete_branch[1] - delete_branch[0] + 1) * SCALE_X)
        .attr('height', seq_height)
        .attr('opacity', 0.3)
        .style('fill', 'black')
        .style('stroke', 'transparent');
    }
    micro_idx++;
  }

  return seq_svgs;
}

function make_schemes(content_selector, Breaks, Strand, Celltype, Type, combined) {
  var dna_length = REF_SEQ[Strand].length;
  var bar_width = dna_length * SCALE_X;

  var labels_width = 500;

  var scale_height = 20;
  var pcis_height = 20;
  var areas_height = 20;
  var seq_height = 20;
  var bars_1_height = 20;

  var microhomologies = MICROHOMOLOGIES.filter(
    x => ((x['hguide'] == Breaks) && (x['strand'] == Strand) && (x['treatment'] == Celltype))
  );
  if (Type != 'all') {
    microhomologies = microhomologies.filter(
      x => x['microhomology_group'] == Type
    );
  }
  microhomologies = microhomologies.sort(
    (x, y) => {
      if (x['microhomology_group'] < y['microhomology_group']) {
        return -1;
      } else if (x['microhomology_group'] > y['microhomology_group']) {
        return 1;
      } else if (x['pos_left'] < y['pos_left']) {
        return -1; 
      } else if (x['pos_left'] > y['pos_left']) {
        return 1;
      } else if (x['pos_right'] < y['pos_right']) {
        return -1;
      } else if (x['pos_right'] > y['pos_right']) {
        return 1;
      } else {
        return 0;
      }
    }
  );

  var scale_y = 0;
  var pcis_y = scale_y + scale_height;
  var areas_y = pcis_y;
  var seq_y = areas_y + areas_height;
  var bars_1_y;
  if (combined) {
    bars_1_y = seq_y + microhomologies.length * seq_height; 
  } else {
    bars_1_y = seq_y + seq_height;
  }
  var scheme_height = bars_1_y + bars_1_height;
  var padding = 10;

  d3.select(content_selector).html("");

  var svgs = combined ?
    init_svgs_combined(content_selector, labels_width, bar_width, scheme_height, padding) :
    init_svgs_separate(content_selector, microhomologies, labels_width, bar_width, scheme_height, padding);

  // The label SVGs
  if (!combined) {
    svgs.append('svg')
      .attr('x', 0)
      .attr('y', 0)
      .attr('width', labels_width)
      .attr('height', scheme_height)
      .style('overflow', 'visible')
      .append('text')
      .attr('width', labels_width)
      .attr('height', scheme_height)
      .attr("dy", scheme_height / 2)
      .style('font-family', 'monospace')
      .text(d => d['microhomology_id']);
  }

  var scale_svgs = svgs.append('svg')
    .attr('x', labels_width)
    .attr('y', scale_y)
    .attr('width', bar_width)
    .attr('height', scale_height);

  var scale = d3.scaleLinear().domain([1, dna_length]).range([SCALE_X / 2, bar_width - SCALE_X / 2]);
  var x_axis = d3.axisBottom().scale(scale).ticks(40);
  scale_svgs.append('g').call(x_axis);

  var pcis_svgs = svgs.append('svg')
    .attr('x', labels_width)
    .attr('y', pcis_y)
    .attr('width', bar_width)
    .attr('height', pcis_height);

  for (var [label, start, end] of PCIS[Strand]) {
    pcis_svgs.append('rect')
      .attr('x', (start - 1) * SCALE_X)
      .attr('y', 0)
      .attr('width', (end - start + 1) * SCALE_X)
      .attr('height', pcis_height)
      .style('fill', 'transparent')
      .style('stroke', PCIS_COLOR)
      .style('stroke-width', 2);

    pcis_svgs.append('text')
      .attr('text-anchor', 'middle')
      .attr('x', (((start + end) / 2) - 1) * SCALE_X)
      .attr('y', 0)
      .attr('dy', '1em')
      .attr('height', pcis_height)
      .text(PCIS_TEXT[label]);
  }

  var areas_svgs = svgs.append('svg')
    .attr('x', labels_width)
    .attr('y', areas_y)
    .attr('width', bar_width)
    .attr('height', areas_height);

  for (var area in AREAS[Breaks][Strand]) {
    if (area == 'intron') {
      continue;
    }

    for (var [start, end] of AREAS[Breaks][Strand][area]) {
      areas_svgs.append('rect')
        .attr('x', (start - 1) * SCALE_X)
        .attr('y', 0)
        .attr('width', (end - start + 1) * SCALE_X)
        .attr('height', areas_height)
        .style('fill', 'transparent')
        .style('stroke', AREAS_COLOR[area])
        .style('stroke-width', 2);

      areas_svgs.append('text')
        .attr('text-anchor', 'middle')
        .attr('x', (((start + end) / 2) - 1) * SCALE_X)
        .attr('y', 0)
        .attr('dy', '1em')
        .attr('height', areas_height)
        .text(AREAS_TEXT[area]);
    }
  }

  if (combined) {
    make_seq_svgs_combined(
      svgs,
      labels_width,
      bar_width,
      seq_y,
      seq_height,
      REF_SEQ[Strand],
      microhomologies,
      DELETE_BRANCH[Celltype] && AREAS[Breaks][Strand]['branch'][0],
    );
  } else {
    make_seq_svgs_separate(
      svgs,
      labels_width,
      bar_width,
      seq_y,
      seq_height,
      REF_SEQ[Strand],
      DELETE_BRANCH[Celltype] && AREAS[Breaks][Strand]['branch'][0],
    );
  }

  var bars_1_svgs = svgs.append('svg')
    .style('overflow', 'visible')
    .attr('x', labels_width)
    .attr('y', bars_1_y)
    .attr('width', bar_width)
    .attr('height', bars_1_height);

  for (var bar of BARS[Strand]) {
    var [bar_label, bar_start, bar_end] = bar;
    var bar_x = (bar_start - 1) * SCALE_X;
    // var bar_y = seq_height;
    var bar_width = (bar_end - bar_start + 1) * SCALE_X;
    var color = BAR_COLOR[bar_label];
    bars_1_svgs.append('rect')
      .attr('x', bar_x)
      .attr('y', 0)
      .attr('width', bar_width)
      .attr('height', bars_1_height)
      .style('fill', color);
    bars_1_svgs.append('text')
      .attr('text-anchor', 'middle')
      .attr('x', bar_x + bar_width / 2)
      .attr('y', 0)
      .attr('dy', '1em')
      .attr('height', bars_1_height)
      .text(BAR_TEXT[bar_label]);
  }
}
