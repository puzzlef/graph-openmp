const fs = require('fs');
const os = require('os');
const path = require('path');

const ROMPTH = /^OMP_NUM_THREADS=(\d+)/m;
const RGRAPH = /^Loading graph .*\/(.*?)\.mtx \.\.\./m;
const RORDER = /^order: (\d+) size: (\d+) \[directed\] \{\}/m;
const RRESLT = /^\{(.+?)ms\} (\w+)/m;
const RRESRV = /^\w+: Reserve space\s*=\s*(.+?) ms/m;
const RADDVT = /^\w+: Add vertices\s*=\s*(.+?) ms/m;
const RRSEDG = /^\w+: Reserve edges\s*=\s*(.+?) ms/m;
const RADDED = /^\w+: Add edges\s*=\s*(.+?) ms/m;
const RUPDAT = /^\w+: Update\s*=\s*(.+?) ms/m;




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
  if (ROMPTH.test(ln)) {
    var [, omp_num_threads] = ROMPTH.exec(ln);
    state.omp_num_threads   = parseFloat(omp_num_threads);
  }
  else if (RGRAPH.test(ln)) {
    var [, graph] = RGRAPH.exec(ln);
    if (!data.has(graph)) data.set(graph, []);
    state.graph = graph;
    state.order = 0;
    state.size  = 0;
    state.reserve_space_time = 0;
    state.add_vertices_time  = 0;
    state.reserve_edges_time = 0;
    state.add_edges_time     = 0;
    state.update_time        = 0;
  }
  else if (RORDER.test(ln)) {
    var [, order, size] = RORDER.exec(ln);
    state.order = parseFloat(order);
    state.size  = parseFloat(size);
  }
  else if (RRESRV.test(ln)) {
    var [, reserve_space_time] = RRESRV.exec(ln);
    state.reserve_space_time   = parseFloat(reserve_space_time);
  }
  else if (RADDVT.test(ln)) {
    var [, add_vertices_time] = RADDVT.exec(ln);
    state.add_vertices_time   = parseFloat(add_vertices_time);
  }
  else if (RRSEDG.test(ln)) {
    var [, reserve_edges_time] = RRSEDG.exec(ln);
    state.reserve_edges_time   = parseFloat(reserve_edges_time);
  }
  else if (RADDED.test(ln)) {
    var [, add_edges_time] = RADDED.exec(ln);
    state.add_edges_time   = parseFloat(add_edges_time);
  }
  else if (RUPDAT.test(ln)) {
    var [, update_time] = RUPDAT.exec(ln);
    state.update_time   = parseFloat(update_time);
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
