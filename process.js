const fs = require('fs');
const os = require('os');
const path = require('path');

const ROMPTH = /^OMP_NUM_THREADS=(\d+)/m;
const RGRAPH = /^Loading graph .*\/(.*?)\.mtx \.\.\./m;
const REXPER = /^Running experiment (\d+) .../m;
const RORDER = /^order: (\d+) size: (\d+) \[directed\] \{\}/m;
const RBATCH = /^Batch fraction: (.+)/m;
const RRESLT = /^\{(.+?)ms(?:; (.+?)ms duplicate)?\} (.+)/m;




// *-FILE
// ------

function readFile(pth) {
  var d = fs.readFileSync(pth, 'utf8');
  return d.replace(/\r?\n/g, '\n');
}

function writeFile(pth, d) {
  d = d.replace(/\r?\n/g, os.EOL);
  fs.writeFileSync(pth, d);
}




// *-CSV
// -----

function writeCsv(pth, rows) {
  var cols = rowKeys(rows[0]);
  var a = cols.join()+'\n';
  for (var r of rows)
    a += [...rowValues(r)].map(v => `"${v}"`).join()+'\n';
  writeFile(pth, a);
}

function rowKeys(row) {
  return Object.keys(row).filter(k => k!=='experiment');
}

function rowValues(row) {
  var a = [];
  for (var k of Object.keys(row))
    if (k!=='experiment') a.push(row[k]);
  return a;
}




// *-LOG
// -----

function readLogLine(ln, data, state) {
  ln = ln.replace(/^\d+-\d+-\d+ \d+:\d+:\d+\s+/, '');
  if (ROMPTH.test(ln)) {
    var [, omp_num_threads] = ROMPTH.exec(ln);
    state.omp_num_threads   = parseFloat(omp_num_threads);
    state.experiment = 0;
    state.graph = '';
    state.order = 0;
    state.size  = 0;
    state.batch_fraction = 0;
    state.time = 0;
    state.duplicate_time = 0;
    state.technique = '';
  }
  else if (REXPER.test(ln)) {
    var [, experiment] = REXPER.exec(ln);
    state.experiment   = parseFloat(experiment);
  }
  else if (RGRAPH.test(ln)) {
    var [, graph] = RGRAPH.exec(ln);
    if (!data.has(graph)) data.set(graph, []);
    state.graph = graph;
  }
  else if (RORDER.test(ln)) {
    var [, order, size] = RORDER.exec(ln);
    state.order = parseFloat(order);
    state.size  = parseFloat(size);
  }
  else if (RBATCH.test(ln)) {
    var [, batch_fraction] = RBATCH.exec(ln);
    state.batch_fraction   = parseFloat(batch_fraction);
  }
  else if (RRESLT.test(ln)) {
    var [, time, duplicate_time, technique] = RRESLT.exec(ln);
    if (technique==='transposeArenaOmp') state.batch_fraction = 0;
    if (technique==='visitCountBfs') {
      var last   = data.get(state.graph).slice(-1)[0];
      technique += last.technique.startsWith('add')? '+' : '-';
    }
    data.get(state.graph).push(Object.assign({}, state, {
      time:           parseFloat(time),
      duplicate_time: parseFloat(duplicate_time || '0'),
      technique,
    }));
  }
  return state;
}

function readLog(pth) {
  var text  = readFile(pth);
  var lines = text.split('\n');
  var data  = new Map();
  var state = {};
  for (var ln of lines)
    state = readLogLine(ln, data, state);
  return data;
}




// PROCESS-*
// ---------

function processCsv(data) {
  var a = [];
  for (var rows of data.values()) {
    for (var r of rows)
      if (r.experiment != 1) a.push(r);
  }
  return a;
}




// MAIN
// ----

function main(cmd, log, out) {
  var data = readLog(log);
  if (path.extname(out)==='') cmd += '-dir';
  switch (cmd) {
    case 'csv':
      var rows = processCsv(data);
      writeCsv(out, rows);
      break;
    case 'csv-dir':
      for (var [graph, rows] of data)
        writeCsv(path.join(out, graph+'.csv'), rows);
      break;
    default:
      console.error(`error: "${cmd}"?`);
      break;
  }
}
main(...process.argv.slice(2));
