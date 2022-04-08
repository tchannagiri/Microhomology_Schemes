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

var WINDOW_NUCLEOTIDES = 10;
var PATTERN_COLOR = '#FFD700';
var WINDOW_COLOR = '#00FF00';
var DELETION_COLOR = 'black';
var HIGHLIGHT_OPACITY = 0.3;

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


function draw_single_seq(
  bar_width,
  ref_seq,
  seq_height,
  seq_svgs,
  microhomology_y,
  microhomology,
  pass, // use multiple passes to make certain elements draw on top
) {
  for (var pos of ['Left', 'Right']) {
    if (pass == 1) {
      // Left or right of match window
      if (pos == 'Left') {
        seq_svgs.append('rect')
          .attr('x', (microhomology[pos] - 1 - WINDOW_NUCLEOTIDES) * SCALE_X)
          .attr('y', microhomology_y)
          .attr('width', WINDOW_NUCLEOTIDES * SCALE_X)
          .attr('height', seq_height)
          .style('fill', WINDOW_COLOR)
          .style('opacity', HIGHLIGHT_OPACITY);
      } else {
        seq_svgs.append('rect')
          .attr('x', (microhomology[pos] - 1 + microhomology['Pattern'].length) * SCALE_X)
          .attr('y', microhomology_y)
          .attr('width', WINDOW_NUCLEOTIDES * SCALE_X)
          .attr('height', seq_height)
          .style('fill', WINDOW_COLOR)
          .style('opacity', HIGHLIGHT_OPACITY);
      }

      // Match nucleotides fill
      seq_svgs.append('rect')
        .attr('x', (microhomology[pos] - 1) * SCALE_X)
        .attr('y', microhomology_y)
        .attr('width', microhomology['Pattern'].length * SCALE_X)
        .attr('height', seq_height)
        .style('fill', PATTERN_COLOR)
        .style('opacity', HIGHLIGHT_OPACITY);
    }

    if (pass == 2) {
      // Match nucleotides outline
      seq_svgs.append('rect')
        .attr('x', (microhomology[pos] - 1) * SCALE_X)
        .attr('y', microhomology_y)
        .attr('width', microhomology['Pattern'].length * SCALE_X)
        .attr('height', seq_height)
        .style('fill', 'transparent')
        .style('stroke', PATTERN_COLOR)
        .style('stroke-width', 2);
    }
  }

  if (pass == 1) {
    // Deletion highlight
    seq_svgs.append('rect')
      .attr('x', (microhomology['Left'] + microhomology['Pattern'].length - 1) * SCALE_X)
      .attr('y', microhomology_y)
      .attr('width', (microhomology['Right'] - microhomology['Left'] - microhomology['Pattern'].length) * SCALE_X)
      .attr('height', seq_height)
      .attr('opacity', HIGHLIGHT_OPACITY)
      .style('fill', DELETION_COLOR);
  }

  if (pass == 2) {
    // Reference sequence text
    seq_svgs.append('text')
      .attr('x', 0)
      .attr('y', microhomology_y)
      .attr("dy", "1em")
      .attr('height', seq_height)
      .attr('textLength', bar_width)
      .style('font-family', 'monospace')
      .style('font-weight', 'bold')
      .text(ref_seq);
  }
}

function make_seq_svgs_separate(
  svgs,
  labels_width,
  bar_width,
  seq_y,
  seq_height,
  ref_seq,
  // delete_branch,
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

  for (var pos of ['Left', 'Right']) {
    seq_svgs.append('rect')
      .attr('x', d => (d[pos] - 1) * SCALE_X)
      .attr('y', 0)
      .attr('width', d => d['Pattern'].length * SCALE_X)
      .attr('height', seq_height)
      .style('fill', 'transparent')
      .style('stroke', PATTERN_COLOR)
      .style('stroke-width', 2);
  }

  seq_svgs.append('rect')
    .attr('x', d => (d['Left'] + d['Pattern'].length - 1) * SCALE_X)
    .attr('y', 0)
    .attr('width', d => (d['Right'] - d['Left'] - d['Pattern'].length) * SCALE_X)
    .attr('height', seq_height)
    .attr('opacity', HIGHLIGHT_OPACITY)
    .style('fill', DELETION_COLOR)
    .style('stroke', 'transparent');

  // if (delete_branch) {
  //   seq_svgs.filter(d => (d['Left'] > delete_branch[1]) || (d['Right'] < delete_branch[0]))
  //     .append('rect')
  //     .attr('x', (delete_branch[0] - 1) * SCALE_X)
  //     .attr('y', 0)
  //     .attr('width', (delete_branch[1] - delete_branch[0] + 1) * SCALE_X)
  //     .attr('height', seq_height)
  //     .attr('opacity', 0.3)
  //     .style('fill', 'black')
  //     .style('stroke', 'transparent');
  // }

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
  // delete_branch,
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

  for (var pass of [1, 2]) {
    for (var micro_idx = 0; micro_idx <  microhomologies.length; micro_idx++) {
      draw_single_seq(
        bar_width,
        ref_seq,
        seq_height,
        seq_svgs,
        micro_idx * seq_height,
        microhomologies[micro_idx],
        pass,
      )
    }
  }
    //   // Left or right of match window
    //   if (pos == 'Left') {
    //     seq_svgs.append('rect')
    //       .attr('x', (micro[pos] - 1 - WINDOW_NUCLEOTIDES) * SCALE_X)
    //       .attr('y', micro_idx * seq_height)
    //       .attr('width', WINDOW_NUCLEOTIDES * SCALE_X)
    //       .attr('height', seq_height)
    //       .style('fill', '#00FF00')
    //       .style('opacity', '0.3');
    //   } else {
    //     seq_svgs.append('rect')
    //       .attr('x', (micro[pos] - 1 + micro['Pattern'].length) * SCALE_X)
    //       .attr('y', micro_idx * seq_height)
    //       .attr('width', WINDOW_NUCLEOTIDES * SCALE_X)
    //       .attr('height', seq_height)
    //       .style('fill', '#00FF00')
    //       .style('opacity', '0.3');
    //   }

    //   // The match nucleotides
    //   // yellow fill
    //   seq_svgs.append('rect')
    //     .attr('x', (micro[pos] - 1) * SCALE_X)
    //     .attr('y', micro_idx * seq_height)
    //     .attr('width', micro['Pattern'].length * SCALE_X)
    //     .attr('height', seq_height)
    //     .style('fill', '#FFFF00')
    //     .style('opacity', '0.3');

    //   // The match nucleotides
    //   // black outline
    //   seq_svgs.append('rect')
    //     .attr('x', (micro[pos] - 1) * SCALE_X)
    //     .attr('y', micro_idx * seq_height)
    //     .attr('width', micro['Pattern'].length * SCALE_X)
    //     .attr('height', seq_height)
    //     .style('fill', 'transparent')
    //     .style('stroke', 'black')
    //     .style('stroke-width', 2);
    // }

    // seq_svgs.append('rect')
    //   .attr('x', (micro['Left'] + micro['Pattern'].length - 1) * SCALE_X)
    //   .attr('y', micro_idx * seq_height)
    //   .attr('width', (micro['Right'] - micro['Left'] - micro['Pattern'].length) * SCALE_X)
    //   .attr('height', seq_height)
    //   .attr('opacity', 0.3)
    //   .style('fill', 'black')
    //   .style('stroke', 'transparent');

    // if (
    //   delete_branch &&
    //   ((micro['Left'] > delete_branch[1]) || (micro['Right'] < delete_branch[0]))
    // ) {
    //   seq_svgs.append('rect')
    //     .attr('x', (delete_branch[0] - 1) * SCALE_X)
    //     .attr('y', micro_idx * seq_height)
    //     .attr('width', (delete_branch[1] - delete_branch[0] + 1) * SCALE_X)
    //     .attr('height', seq_height)
    //     .attr('opacity', 0.3)
    //     .style('fill', 'black')
    //     .style('stroke', 'transparent');
    // }
    // micro_idx++;

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
      .text(microhomologies[i]['Name']);

    // seq_svgs.append('text')
    //   .attr('x', 0)
    //   .attr('y', i * seq_height)
    //   .attr("dy", "1em")
    //   .attr('height', seq_height)
    //   .attr('textLength', bar_width)
    //   .style('font-family', 'monospace')
    //   .text(ref_seq);
  }

  return seq_svgs;
}

function make_schemes(content_selector, Breaks, Strand, Celltype, Type, combined) {
  var dna_length = REF_SEQ[Celltype][Strand].length;
  var bar_width = dna_length * SCALE_X;

  var labels_width = 100;

  var pad = 3;
  var scale_height = 20;
  var pcis_height = 20;
  var areas_height = 20;
  var seq_height = 20;
  var bars_1_height = 20;

  var microhomologies = MICROHOMOLOGIES.filter(
    x => ((x['Breaks'] == Breaks) && (x['Strand'] == Strand) && (x['Celltype'] == Celltype))
  );
  if (Type != 'all') {
    microhomologies = microhomologies.filter(
      x => x['Type'] == Type
    );
  }
  // microhomologies = microhomologies.sort(
  //   (x, y) => {
  //     if (x['microhomology_group'] < y['microhomology_group']) {
  //       return -1;
  //     } else if (x['microhomology_group'] > y['microhomology_group']) {
  //       return 1;
  //     } else if (x['pos_left'] < y['pos_left']) {
  //       return -1; 
  //     } else if (x['pos_left'] > y['pos_left']) {
  //       return 1;
  //     } else if (x['pos_right'] < y['pos_right']) {
  //       return -1;
  //     } else if (x['pos_right'] > y['pos_right']) {
  //       return 1;
  //     } else {
  //       return 0;
  //     }
  //   }
  // );

  var scale_y = 0;
  var pcis_y = scale_y + scale_height + pad;
  var areas_y = pcis_y;
  var seq_y = areas_y + areas_height + pad;
  var bars_1_y;
  if (combined) {
    bars_1_y = seq_y + microhomologies.length * seq_height + pad; 
  } else {
    bars_1_y = seq_y + seq_height + pad;
  }
  var scheme_height = bars_1_y + bars_1_height + pad;
  var seq_pad = 10;

  d3.select(content_selector).html("");

  var svgs = combined ?
    init_svgs_combined(content_selector, labels_width, bar_width, scheme_height, seq_pad) :
    init_svgs_separate(content_selector, microhomologies, labels_width, bar_width, scheme_height, seq_pad);

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
      .text(d => d['Name']);
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

  for (var pcis of PCIS[Celltype][Strand]) {
    pcis_svgs.append('rect')
      .attr('x', (pcis['start'] - 1) * SCALE_X)
      .attr('y', 0)
      .attr('width', (pcis['end'] - pcis['start'] + 1) * SCALE_X)
      .attr('height', pcis_height)
      .style('fill', 'transparent')
      .style('stroke', PCIS_COLOR)
      .style('stroke-width', 2);

    pcis_svgs.append('text')
      .attr('text-anchor', 'middle')
      .attr('x', (((pcis['start'] + pcis['end']) / 2) - 1) * SCALE_X)
      .attr('y', 0)
      .attr('dy', '1em')
      .attr('height', pcis_height)
      .text(PCIS_TEXT[pcis['name']]);
  }

  var areas_svgs = svgs.append('svg')
    .attr('x', labels_width)
    .attr('y', areas_y)
    .attr('width', bar_width)
    .attr('height', areas_height);

  for (var area of AREAS[Celltype][Breaks][Strand]) {
    areas_svgs.append('rect')
      .attr('x', (area['start'] - 1) * SCALE_X)
      .attr('y', 0)
      .attr('width', (area['end'] - area['start'] + 1) * SCALE_X)
      .attr('height', areas_height)
      .style('fill', 'transparent')
      .style('stroke', AREAS_COLOR[area['name']])
      .style('stroke-width', 2);

    areas_svgs.append('text')
      .attr('text-anchor', 'middle')
      .attr('x', (((area['start'] + area['end']) / 2) - 1) * SCALE_X)
      .attr('y', 0)
      .attr('dy', '1em')
      .attr('height', areas_height)
      .text(AREAS_TEXT[area['name']]);
  }

  if (combined) {
    make_seq_svgs_combined(
      svgs,
      labels_width,
      bar_width,
      seq_y,
      seq_height,
      REF_SEQ[Celltype][Strand],
      microhomologies,
      // DELETE_BRANCH[Celltype],
    );
  } else {
    make_seq_svgs_separate(
      svgs,
      labels_width,
      bar_width,
      seq_y,
      seq_height,
      REF_SEQ[Celltype][Strand],
      DELETE_BRANCH[Celltype],
    );
  }

  var bars_1_svgs = svgs.append('svg')
    .style('overflow', 'visible')
    .attr('x', labels_width)
    .attr('y', bars_1_y)
    .attr('width', bar_width)
    .attr('height', bars_1_height);

  for (var bar of BARS[Celltype][Strand]) {
    var bar_x = (bar['start'] - 1) * SCALE_X;
    var bar_width = (bar['end'] - bar['start'] + 1) * SCALE_X;
    var color = BAR_COLOR[bar['name']];
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
      .text(BAR_TEXT[bar['name']]);
  }
}
