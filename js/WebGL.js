let canvas, ctx;

// get the canvas element using the DOM
canvas = document.getElementById('mycanvas');

let gl = canvas.getContext('webgl');

gl.clearColor(0.9, 0.9, 0.8, 1);
gl.clear(gl.COLOR_BUFFER_BIT);