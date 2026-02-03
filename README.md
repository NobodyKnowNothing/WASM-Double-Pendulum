# WASM Double Pendulum Simulator

A high-performance double pendulum simulation utilizing **WebAssembly (WASM)** for numerical integration. This project demonstrates the power of WASM by running complex physics simulations directly in the browser with near-native performance.

The simulator compares three different integration methods side-by-side:
1.  **Naive Runge-Kutta 4 (RK4)**: Standard fixed-step integration.
2.  **Adaptive RK4**: integration with adaptive step sizing for error control.
3.  **Boost.Numeric.ODEInt (RKD5)**: Using the simpler Runge-Kutta-Dopri5 stepper from the Boost C++ libraries.

## Features

*   **Real-time Visualization**: Double pendulum animation using P5.js.
*   **Live Data Plotting**: Real-time graphing of Lyapunov Exponents and Total Energy using Plotly.js.
*   **Chaos Detection**: Calculates Lyapunov exponents to quantify the chaotic nature of the system.
*   **Interactive Controls**: Modify mass, lengths, gravity, inititial conditions, and simulation parameters on the fly.
*   **Comparison Mode**: Visual separation of different integration algorithms to observe divergence over time.

## Project Structure

WASM-Double-Pendulum/
├── src/
│   └── cpp/                # C++ Source Code
│       ├── main_web.cpp
│       └── simulation.hpp
├── public/                 # Web Application (Serve this folder)
│   ├── index.html          # Entry point (updated paths)
│   ├── js/
│   │   └── simAnimator.js  # Visualization logic
│   └── wasm/
│       ├── simulation.js   # Emscripten glue code
│       └── simulation.wasm # Compiled WebAssembly binary
└── README.md               # Project documentation

## Getting Started

### Prerequisites

*   A modern web browser with WebAssembly support.
*   A local web server (to bypass CORS restrictions when loading WASM modules).

### Running the Simulator

1.  Clone the repository.
2.  Navigate to the project root.
3.  Serve the `public` directory.
    *   **Python**: `cd public && python -m http.server`
    *   **Node (http-server)**: `npx http-server public`
    *   **VS Code**: Open `public/index.html` and use the "Live Server" extension.
4.  Open `localhost:8000` (or the port provided by your server) in your browser.

## Building from Source

To modify the C++ simulation code and recompile the WASM module:

1.  Install the [Emscripten SDK (emsdk)](https://emscripten.org/docs/getting_started/downloads.html).
2.  Activate the emsdk environment.
3.  Run the compile command (example):

```sh
emcc src/cpp/main_web.cpp -I src/cpp -o public/wasm/simulation.js \
    -O3 \
    -s EXPORTED_FUNCTIONS="['_simulation', '_main', '_malloc', '_free']" \
    -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap']" \
    -s ALLOW_MEMORY_GROWTH=1
```

*Note: Ensure you include Boost headers in your include path if they are not in the standard emscripten include path, typically using `-I path/to/boost`.*

## Technologies Used

*   **C++**: Core physics and integration logic.
*   **Emscripten**: Compiling C++ to WebAssembly.
*   **Boost Libraries**: Advanced numerical integration (`boost::numeric::odeint`).
*   **JavaScript**: Canvas rendering and UI logic.
*   **P5.js**: Animation and drawing.
*   **Plotly.js**: Data visualization.
