const fs = require('fs');
const os = require('os');
const path = require('path');

const RGRAPH = /^Loading graph .*\/(.*?)\.mtx \.\.\./m;
const ROMPTH = /^OMP_NUM_THREADS=(\d+)/m;
const RPARTS = /^readMtxOmpW\(\): vertices=([\d.]+)ms, read=([\d.]+)ms, parse=([\d.]+)ms, edges=([\d.]+)ms, update=([\d.]+)ms/;
const RORDER = /^order: (\d+) size: (\d+) (?:\[\w+\] )?\{\}/m;
const RRESLT = /^\[(\S+?) ms\] \w+/m;




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
  var cols = Object.keys(rows[0]);
  var a = cols.join()+'\n';
  for (var r of rows)
    a += [...Object.values(r)].map(v => `"${v}"`).join()+'\n';
  writeFile(pth, a);
}




// *-LOG
// -----

function readLogLine(ln, data, state) {
  ln = ln.replace(/^\d+-\d+-\d+ \d+:\d+:\d+\s+/, '');
  if (RGRAPH.test(ln)) {
    var [, graph] = RGRAPH.exec(ln);
    if (!data.has(graph)) data.set(graph, []);
    state = {graph};
  }
  else if (ROMPTH.test(ln)) {
    var [, omp_num_threads] = ROMPTH.exec(ln);
    state.omp_num_threads   = parseFloat(omp_num_threads);
    state.order = 0;
    state.size  = 0;
    state.vertices_time = 0;
    state.read_time     = 0;
    state.parse_time    = 0;
    state.edges_time    = 0;
    state.update_time   = 0;
  }
  else if (RPARTS.test(ln)) {
    var [, vertices_time, read_time, parse_time, edges_time, update_time] = RPARTS.exec(ln);
    state.vertices_time = parseFloat(vertices_time);
    state.read_time     = parseFloat(read_time);
    state.parse_time    = parseFloat(parse_time);
    state.edges_time    = parseFloat(edges_time);
    state.update_time   = parseFloat(update_time);
  }
  else if (RORDER.test(ln)) {
    var [, order, size] = RORDER.exec(ln);
    state.order = parseFloat(order);
    state.size  = parseFloat(size);
  }
  else if (RRESLT.test(ln)) {
    var [, total_time] = RRESLT.exec(ln);
    data.get(state.graph).push(Object.assign({}, state, {
      total_time: parseFloat(total_time),
    }));
  }
  return state;
}

function readLog(pth) {
  var text  = readFile(pth);
  var lines = text.split('\n');
  var data  = new Map();
  var state = null;
  for (var ln of lines)
    state = readLogLine(ln, data, state);
  return data;
}




// PROCESS-*
// ---------

function processCsv(data) {
  var a = [];
  for (var rows of data.values())
    a.push(...rows);
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
