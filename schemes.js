var BAR_COLOR = {
  'dsred': 'red',
  'intron': 'pink',
};

var BAR_TEXT = {
  'dsred': 'DS Red',
  'intron': 'Intron',
};

var AREAS_COLOR = {
  'exon1': '#4471C4',
  'exon2': '#4471C4',
  'intron': '#00AF50',
  'branch': 'gray',
};

var AREAS_TEXT = {
  'exon1': 'Exon',
  'exon2': 'Exon',
  'intron': 'Intron',
  'branch': 'Intron (branch site)',
};

var PCIS_TEXT = {
  // 'r5': 'DsRed.pCis.R5',
  // 'f6': 'DsRed.pCis.F6',
  'r5': 'Primer',
  'f6': 'Primer',
};

var PCIS_COLOR = 'purple';

var DELETE_BRANCH = {
  'wt': false,
  'db': true,
};

var SCALE_X = 10;
var SCALE_Y = 20;

var WINDOW_NUCLEOTIDES = 10;
var PATTERN_COLOR = '#FFC000';
var WINDOW_COLOR = '#00FF00';
var DELETION_COLOR = 'black';
var HIGHLIGHT_OPACITY = 0.3;

var LABELS = {
  'hg39': 'sgRNA A',
  'hg42': 'sgRNA B',
  '2dsb': 'sgRNA A & B',
  'Forward': 'Forward strand',
  'Reverse': 'Reverse strand',
  'wt': 'WT',
  'db': 'BranchΔ',
}

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
    }

    if (pass == 2) {
      // Match nucleotides fille/outline
      seq_svgs.append('rect')
        .attr('x', (microhomology[pos] - 1) * SCALE_X)
        .attr('y', microhomology_y)
        .attr('width', microhomology['Pattern'].length * SCALE_X)
        .attr('height', seq_height)
        .style('fill', PATTERN_COLOR)
        .style('fill-opacity', HIGHLIGHT_OPACITY)
        .style('stroke', PATTERN_COLOR)
        .style('stroke-opacity', 1)
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
  scheme_svgs,
  labels_width,
  content_width,
  seq_y,
  seq_height,
  ref_seq,
) {
  var seq_svgs = scheme_svgs.append('g')
    .style('overflow', 'visible')
    .attr('transform', `translate(${labels_width}, ${seq_y})`)
    .attr('width', content_width)
    .attr('height', seq_height);
    
    for (var pos of ['Left', 'Right']) {
      // Window highlights
      if (pos == 'Left') {
        seq_svgs.append('rect')
          .attr('x', d => (d[pos] - 1 - WINDOW_NUCLEOTIDES) * SCALE_X)
          .attr('y', 0)
          .attr('width', WINDOW_NUCLEOTIDES * SCALE_X)
          .attr('height', seq_height)
          .style('fill', WINDOW_COLOR)
          .style('opacity', HIGHLIGHT_OPACITY);
      } else {
        seq_svgs.append('rect')
          .attr('x', d => (d[pos] + d['Pattern'].length - 1) * SCALE_X)
          .attr('y', 0)
          .attr('width', WINDOW_NUCLEOTIDES * SCALE_X)
          .attr('height', seq_height)
          .style('fill', WINDOW_COLOR)
          .style('opacity', HIGHLIGHT_OPACITY);
      }

      // Pattern highlights
      seq_svgs.append('rect')
        .attr('x', d => (d[pos] - 1) * SCALE_X)
        .attr('y', 0)
        .attr('width', d => d['Pattern'].length * SCALE_X)
        .attr('height', seq_height)
        .style('fill', PATTERN_COLOR)
        .style('opacity', HIGHLIGHT_OPACITY);
    }

  // Pattern outline
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

  // Deletion highlight
  seq_svgs.append('rect')
    .attr('x', d => (d['Left'] + d['Pattern'].length - 1) * SCALE_X)
    .attr('y', 0)
    .attr('width', d => (d['Right'] - d['Left'] - d['Pattern'].length) * SCALE_X)
    .attr('height', seq_height)
    .attr('opacity', HIGHLIGHT_OPACITY)
    .style('fill', DELETION_COLOR)
    .style('stroke', 'transparent');
  
  // Reference sequence text
  seq_svgs.append('text')
    .attr('x', 0)
    .attr('y', seq_height / 2)
    .attr('dy', '.5ex')
    .attr('textLength', content_width)
    .style('font-family', 'monospace')
    .style('font-weight', 'bold')
    .text(ref_seq);
    
  return seq_svgs;
}

function make_seq_svgs_combined(
  scheme_svg,
  labels_width,
  bar_width,
  seq_y,
  seq_height,
  ref_seq,
  microhomologies,
  cut_pos,
) {
  var label_svgs = scheme_svg.append('g')
    .style('overflow', 'visible')
    .attr('transform', `translate(0, ${seq_y})`)
    .attr('width', labels_width)
    .attr('height', seq_height * microhomologies.length);

  var seq_svgs = scheme_svg.append('g')
    .style('overflow', 'visible')
    .attr('transform', `translate(${labels_width}, ${seq_y})`)
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

  // The microhomology names
  for (var micro_idx = 0; micro_idx < microhomologies.length; micro_idx++) {
    label_svgs.append('g')
        .attr('transform', `translate(0, ${micro_idx * seq_height})`)
        .attr('width', labels_width)
        .attr('height', seq_height)
        .style('overflow', 'visible')
      .append('text')
        .attr('x', 0)
        .attr('y', seq_height / 2)
        .attr('dy', '.5ex')
        .style('font-family', 'monospace')
        .style('font-weight', microhomologies[micro_idx]['Bold'] == 1 ? 'bold' : 'normal')
        .text(microhomologies[micro_idx]['Name'] + (microhomologies[micro_idx]['Bold'] == 1 ? ' *' : ''));
  }

  return seq_svgs;
}

function get_title(Breaks, Strand, Celltype, Type) {
  var title = [
    LABELS[Breaks], 
    LABELS[Strand], 
    LABELS[Celltype], 
    LABELS[Type], 
  ];
  title = title.filter(x => x != null);
  title = title.join(', ')
  return 'Microhomology schemes: ' + title;
}

function make_schemes(content_selector, Breaks, Strand, Celltype, Type, combined) {
  var dna_length = REF_SEQ[Celltype][Strand].length;
  var scheme_width = dna_length * SCALE_X;
  var labels_width = 100;
  var total_width = labels_width + scheme_width;

  var title_font_size = 30;

  var areas_stroke_width = 3;

  var pad = 3;
  var title_height = 50;
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
  var bars_1_y = scale_y + scale_height + pad; 
  var pcis_y = bars_1_y + bars_1_height + pad;
  var areas_y = pcis_y;
  var seq_y =  pcis_y + pcis_height + pad; 

  var scheme_height;
  var total_scheme_height;
  if (combined) {
    scheme_height = null;
    total_scheme_height = seq_y + microhomologies.length * seq_height; 
  } else {
    scheme_height = seq_y + seq_height + pad;
    total_scheme_height = scheme_height * microhomologies.length;
  }

  var title_y = 0;
  var scheme_y = title_y + title_height + pad;
  var main_height = scheme_y + total_scheme_height;

  d3.select(content_selector).html('');
  var main_svg = d3.select(content_selector).append('svg')
    .attr('width', total_width)
    .attr('height', main_height);

  // The title
  main_svg.append('g')
      .attr('width', total_width)
      .attr('height', title_height)
    .append('text')
      .attr('x', total_width / 2)
      .attr('y', title_height / 2)
      .style('text-anchor', 'middle')
      .style('dominant-baseline', 'middle')
      .style('overflow', 'visible')
      .style('font-size', title_font_size)
      .text(get_title(Breaks, Strand, Celltype, Type));

  var scheme_svg_outer = main_svg.append('g')
    .attr('transform', `translate(0, ${scheme_y})`)
    .attr('width', total_width)
    .attr('height', total_scheme_height);

  var scheme_svg;
  if (combined) {
    scheme_svg = scheme_svg_outer;
  } else {
    scheme_svg = scheme_svg_outer
      .selectAll('g')
      .data(microhomologies)
      .enter()
      .append('g')
        .attr('transform', (d, i) => `translate(0, ${i * scheme_height})`)
        .attr('width', total_width)
        .attr('height', scheme_height);
  }

  // The label SVGs
  if (!combined) {
    scheme_svg.append('g')
        .attr('width', labels_width)
        .attr('height', scheme_height)
        .style('overflow', 'visible')
      .append('text')
        .attr('x', 0)
        .attr('y', scheme_height / 2)
        .attr('dy', '.5ex')
        .style('font-family', 'monospace')
        .style('font-weight', d => d['Bold'] == 1 ? 'bold' : 'light')
        .text(d => d['Name'] + (d['Bold'] == 1 ? ' *' : ''));
  }

  // The scale bar at the top
  var scale_svgs = scheme_svg.append('g')
    .attr('transform', `translate(${labels_width}, ${scale_y})`)
    .attr('width', scheme_width)
    .attr('height', scale_height);

  var scale = d3.scaleLinear().domain([1, dna_length]).range([SCALE_X / 2, scheme_width - SCALE_X / 2]);
  var x_axis = d3.axisBottom().scale(scale).ticks(40);
  scale_svgs.append('g').call(x_axis);

  // The bars for DsRed and Intron
  var bars_1_svgs = scheme_svg.append('g')
    .style('overflow', 'visible')
    .attr('transform', `translate(${labels_width}, ${bars_1_y})`)
    .attr('width', scheme_width)
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
      .attr('y', bars_1_height / 2)
      .attr('dy', '.5ex')
      .text(BAR_TEXT[bar['name']]);
  }

  // The Primer windows
  var pcis_svgs = scheme_svg.append('g')
    .attr('transform', `translate(${labels_width}, ${pcis_y})`)
    .attr('width', scheme_width)
    .attr('height', pcis_height);

  for (var pcis of PCIS[Celltype][Strand]) {
    pcis_svgs.append('rect')
      .attr('x', (pcis['start'] - 1) * SCALE_X)
      .attr('y', 0)
      .attr('width', (pcis['end'] - pcis['start'] + 1) * SCALE_X)
      .attr('height', pcis_height)
      .style('fill', 'rgba(0, 0, 0, 0)')
      .style('fill-opacity', 0)
      .style('stroke', PCIS_COLOR)
      .style('stroke-opacity', 1)
      .style('stroke-width', areas_stroke_width);

    pcis_svgs.append('text')
      .attr('text-anchor', 'middle')
      .attr('x', (((pcis['start'] + pcis['end']) / 2) - 1) * SCALE_X)
      .attr('y', pcis_height / 2)
      .attr('dy', '.5ex')
      .text(PCIS_TEXT[pcis['name']]);
  }

  // The Intron/Exon/Branch windows
  var areas_svgs = scheme_svg.append('g')
    .attr('transform', `translate(${labels_width}, ${areas_y})`)
    .attr('width', scheme_width)
    .attr('height', areas_height);

  for (var area of AREAS[Celltype][Breaks][Strand]) {
    areas_svgs.append('rect')
      .attr('x', (area['start'] - 1) * SCALE_X)
      .attr('y', 0)
      .attr('width', (area['end'] - area['start'] + 1) * SCALE_X)
      .attr('height', areas_height)
      .style('fill', 'rgba(0, 0, 0, 0)')
      .style('fill-opacity', 0)
      .style('stroke', AREAS_COLOR[area['name']])
      .style('stroke-opacity', 1)
      .style('stroke-width', areas_stroke_width);

    areas_svgs.append('text')
      .attr('text-anchor', 'middle')
      .attr('x', (((area['start'] + area['end']) / 2) - 1) * SCALE_X)
      .attr('y', areas_height / 2)
      .attr('dy', '.5ex')
      .text(AREAS_TEXT[area['name']]);
  }

  // The actual microhomologies
  if (combined) {
    make_seq_svgs_combined(
      scheme_svg,
      labels_width,
      scheme_width,
      seq_y,
      seq_height,
      REF_SEQ[Celltype][Strand],
      microhomologies,
    );
  } else {
    make_seq_svgs_separate(
      scheme_svg,
      labels_width,
      scheme_width,
      seq_y,
      seq_height,
      REF_SEQ[Celltype][Strand],
    );
  }
}
