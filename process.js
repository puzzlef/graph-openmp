const fs = require('fs');
const os = require('os');
const path = require('path');

const ROMPTH = /^OMP_NUM_THREADS=(\d+)/m;
const RGRAPH = /^Loading graph .*\/(.*?)\.mtx \.\.\./m;
const RALLOC = /^readMtxFormatToCsrOmpW: \{(.+?)ms\} Allocate memory/m;
const RREADE = /^readMtxFormatToCsrOmpW: \{(.+?)ms\} Read Edgelist/m;
const RCONVC = /^readMtxFormatToCsrOmpW: \{(.+?)ms\} Convert to CSR/m;
const RFREEM = /^readMtxFormatToCsrOmpW: \{(.+?)ms\} Free memory/m;
const RORDER = /^order: (\d+) size: (\d+) \[directed\] \{\}/m;
const RRESLT = /^\{(.+?)ms\} (.+)/m;




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
  return Object.keys(row);
}

function rowValues(row) {
  var a = [];
  for (var k of Object.keys(row))
    a.push(row[k]);
  return a;
}




// *-LOG
// -----

function readLogLine(ln, data, state) {
  ln = ln.replace(/^\d+-\d+-\d+ \d+:\d+:\d+\s+/, '');
  if (ROMPTH.test(ln)) {
    var [, omp_num_threads] = ROMPTH.exec(ln);
    state.omp_num_threads   = parseFloat(omp_num_threads);
    state.graph = '';
    state.order = 0;
    state.size  = 0;
    state.allocation_time = 0;
    state.read_time       = 0;
    state.convert_time    = 0;
    state.free_time       = 0;
    state.time      = 0;
    state.technique = '';
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
  else if (RALLOC.test(ln)) {
    var [, time] = RALLOC.exec(ln);
    state.allocation_time = parseFloat(time);
  }
  else if (RREADE.test(ln)) {
    var [, time] = RREADE.exec(ln);
    state.read_time = parseFloat(time);
  }
  else if (RCONVC.test(ln)) {
    var [, time] = RCONVC.exec(ln);
    state.convert_time = parseFloat(time);
  }
  else if (RFREEM.test(ln)) {
    var [, time] = RFREEM.exec(ln);
    state.free_time = parseFloat(time);
  }
  else if (RRESLT.test(ln)) {
    var [, time, technique] = RRESLT.exec(ln);
    data.get(state.graph).push(Object.assign({}, state, {
      time: parseFloat(time),
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
      a.push(r);
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
