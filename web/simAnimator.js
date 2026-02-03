function memAllocateWasm(params) { // Converts an array into a pointer accessible by wasm
    if (!Module._malloc) { // Checks if module is instantiated 
        console.error("Malloc not yet defined.");
        return null;
    }
    const nDataBytes = params.length * Float64Array.BYTES_PER_ELEMENT; // Defined amount of bytes required

    const dataPtr = Module._malloc(nDataBytes); // Allocates those bytes into WASM

    let buffer = Module.HEAPU8.buffer; // Defines buffer, basically the entire block of ram, Dont ever store as a global variable as memory can be resized for optimization and this buffer object becomes wrong 

    const dataHeap = new Float64Array(buffer, dataPtr, params.length); // Specifies part of ram to access, makes a 'window' into the ram

    dataHeap.set(params); // Sets the window to params

    return dataPtr; // Returns pointer to the window
}

function simWrapper(params) { // Runs simulation with params data
    const paramsPtr = memAllocateWasm(params); // Allocates wasm memory and returns pointer

    console.log("Emscripten Module Object:", Module);

    let resultPtr = Module._simulation(paramsPtr); // Runs simulation with params data from the pointer

    const size = new Float64Array(Module.HEAPU8.buffer, resultPtr, 1)[0]; // Creates a ram window of size 1 at the output data pointer to retreive the size number specified in WASM

    console.log("Total Elements: ", size);

    const resultData = new Float64Array(Module.HEAPU8.buffer, resultPtr, size); // Creates a ram window according to the size previously retreived

    Module._free(paramsPtr); // Frees params pointer 

    return resultData; // t, th1, w1, th2, w2, l1, l2, l3, l4
}

// Global variables for animation

let naiverk4 = null;
let rk4 = null;
let boost = null;
let divider1 = 0;
let divider2 = 0;
let divider3 = 0;
let fps = 24;
let loaded = false;
let l1 = 1;
let l2 = 1;
let th1 = 0;
let th2 = 0;

let x1 = 200;
let y1 = 140;

let x2 = 0;
let y2 = 0;

let x3 = 0;
let y3 = 0;

let radOffset = 0

let index = 0;
let plotCounter = 0;

let stats = {
    naive: { sumLog: 0, time: 0 },
    rk4: { sumLog: 0, time: 0 },
    boost: { sumLog: 0, time: 0 }
};

// Small helper to safely calculate log (prevents errors if value is 0 or negative)
function getSafeLog(val) {
    let absoluteVal = Math.abs(val);
    return Math.log10(absoluteVal + 1e-15); // Use 1e-15 to avoid log(0) = -Infinity
}

function runSimulation() {

    index = 0;
    plotCounter = 0;
    initGraph();
    stats = {
        naive: { sumLog: 0.01, count: 0.01 },
        rk4: { sumLog: 0.01, count: 0.01 },
        boost: { sumLog: 0.01, count: 0.01 }
    };
    if (!loaded) return false;

    // Retreives info from html and passes to global variables required for drawing JS animation instead of just passing to WASM
    l1 = parseFloat(document.getElementById("length-1").value);
    l2 = parseFloat(document.getElementById("length-2").value);
    fps = parseFloat(document.getElementById("fps").value);
    th1 = parseFloat(document.getElementById("init-state-1").value);
    th2 = parseFloat(document.getElementById("init-state-3").value);


    /* 
    Defines input parameter array for WASM
    params [m1, m2, l1, l2, t_start, t_end, steps, g, state[0], state[1], state[2], state[3], fps]
    */
    const params = [
        parseFloat(document.getElementById("mass-1").value),
        parseFloat(document.getElementById("mass-2").value),
        l1,
        l2,
        parseFloat(document.getElementById("t-start").value),
        parseFloat(document.getElementById("t-end").value),
        parseFloat(document.getElementById("steps").value),
        parseFloat(document.getElementById("gravity").value),
        th1,
        parseFloat(document.getElementById("init-state-2").value),
        th2,
        parseFloat(document.getElementById("init-state-4").value),
        fps
    ];

    resultData = simWrapper(params); // Actual WASM Call, passes params array to WASM by passing params array pointer, WASM returns the result data by pointer and this wrapper turns it back into an array 

    // Checks if there are enough steps per second to draw at the specified fps, is not a strict requirement but will cause issues if not met
    if ((params[6] / (params[5] - params[4])) < (params[12])) console.warn("Less steps per second then frames per second, will cause issues.");

    // Splits data into 3 arrays, for each method of integration
    divider3 = resultData[0]; // Defines total length of data, and end index of boost integration
    divider1 = resultData[1]; // Defines end index of naive rk4 integration
    divider2 = resultData[2]; // Defines end index of rk4 integration

    naiverk4 = resultData.slice(3, divider1); // Defines data for naive rk4 integration
    rk4 = resultData.slice(divider1, divider2); // Defines data for rk4 integration
    boost = resultData.slice(divider2, divider3); // Defines data for boost integration

    console.log(naiverk4);
    console.log(rk4);
    console.log(boost);
}

function setup() {
    let canvas = createCanvas(400, 400);

    canvas.parent('canvas-container');

    frameRate(fps);

    radOffset = PI / 2;
    th1 = 0;
    th2 = 0;

    x2 = x1 + round(40 * l1 * cos(th1 + radOffset));
    y2 = y1 + round(40 * l1 * sin(th1 + radOffset));

    x3 = x2 + round(40 * l2 * cos(th2 + radOffset));
    y3 = y2 + round(40 * l2 * sin(th2 + radOffset));

    initGraph();
}

function drawDoubleP(th1, th2, colour) {
    stroke(colour);

    // w1 = random(-1, 1) * 0.1;
    // w2 = random(-1, 1) * 0.1;

    th1 += radOffset;
    th2 += radOffset;

    x2 = x1 + round(40 * l1 * cos(th1));
    y2 = y1 + round(40 * l1 * sin(th1));

    x3 = x2 + round(40 * l2 * cos(th2));
    y3 = y2 + round(40 * l2 * sin(th2));

    line(x1, y1, x2, y2);
    line(x2, y2, x3, y3);
}


indexJump = 10;
function draw() {
    background(220);
    if (naiverk4 == null) return; // || divider == 0

    const naiverk4Line = naiverk4.slice(index, index + indexJump);
    const rk4Line = rk4.slice(index, index + indexJump);
    const boostLine = boost.slice(index, index + indexJump);

    drawDoubleP(naiverk4Line[1], naiverk4Line[3], 'red');
    drawDoubleP(rk4Line[1], rk4Line[3], 'blue');
    drawDoubleP(boostLine[1], boostLine[3], 'green');

    index += indexJump;
    if (index > (divider1 - indexJump)) { index = 0; }

    stats.naive.sumLog = naiverk4Line[5];
    stats.naive.time = naiverk4Line[0];

    stats.rk4.sumLog = rk4Line[5];
    stats.rk4.time = rk4Line[0];

    stats.boost.sumLog = boostLine[5];
    stats.boost.time = boostLine[0];

    plotCounter++;
    if (plotCounter % 1 === 0 && plotCounter < min(divider1 / indexJump, divider2 / indexJump - divider1 / indexJump, divider3 / indexJump - divider2 / indexJump)) {
        Plotly.extendTraces('chart-container1', {
            y: [[stats.naive.sumLog / stats.naive.time], [stats.rk4.sumLog / stats.rk4.time], [stats.boost.sumLog / stats.boost.time]],
            x: [[stats.naive.time], [stats.rk4.time], [stats.boost.time]]
        }, [0, 1, 2], 300);

        Plotly.extendTraces('chart-container2', {
            y: [[naiverk4Line[9]], [rk4Line[9]], [boostLine[9]]],
            x: [[stats.naive.time], [stats.rk4.time], [stats.boost.time]]
        }, [0, 1, 2], 300);
    }
}

function initGraph() {
    initGraph1();
    initGraph2();
}

function initGraph1() {
    let trace1 = {
        x: [],
        y: [],
        mode: 'lines',
        line: { color: 'red', width: 2 },
        name: 'Naive RK4'
    };

    let trace2 = {
        x: [],
        y: [],
        mode: 'lines',
        line: { color: 'blue', width: 2 },
        name: 'RK4'
    };

    let trace3 = {
        x: [],
        y: [],
        mode: 'lines',
        line: { color: 'green', width: 2 },
        name: 'Boost'
    };

    let layout = {
        title: 'Largest Lyapunov Exponent (LLE) vs Time',
        xaxis: { title: 'Time (s)' },
        yaxis: {
            title: 'LLE (λ)',
            autorange: true // Let it scale automatically at first
        },
        // IMPORTANT: staticPlot: true or following makes it faster
        hovermode: false,
        margin: { t: 40, b: 40, l: 50, r: 20 }
    };

    Plotly.newPlot('chart-container1', [trace1, trace2, trace3], layout, { responsive: true, displayModeBar: false });
}

function initGraph2() {
    let trace1 = {
        x: [],
        y: [],
        mode: 'lines',
        line: { color: 'red', width: 2 },
        name: 'Naive RK4'
    };

    let trace2 = {
        x: [],
        y: [],
        mode: 'lines',
        line: { color: 'blue', width: 2 },
        name: 'RK4'
    };

    let trace3 = {
        x: [],
        y: [],
        mode: 'lines',
        line: { color: 'green', width: 2 },
        name: 'Boost'
    };

    let layout = {
        title: 'Total Energy (E) vs Time (s)',
        xaxis: { title: 'Time (s)' },
        yaxis: {
            title: 'Energy (J)',
            autorange: true // Let it scale automatically at first
        },
        // IMPORTANT: staticPlot: true or following makes it faster
        hovermode: false,
        margin: { t: 40, b: 40, l: 50, r: 20 }
    };

    Plotly.newPlot('chart-container2', [trace1, trace2, trace3], layout, { responsive: true, displayModeBar: false });
}

function checkModuleStatus() {
    if (typeof Module !== 'undefined' && Module.calledRun) {
        console.log("Module already loaded");
        loaded = true;
        runSimulation();
    } else {
        // Fallback: overwrite the callback
        Module.onRuntimeInitialized = () => {
            console.log("Module initialized via callback");
            loaded = true;
            runSimulation();
        };
    }
}
window.addEventListener('DOMContentLoaded', () => {
    checkModuleStatus();
});
