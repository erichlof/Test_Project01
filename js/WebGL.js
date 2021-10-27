let canvas; 
let ctx;
let gl;
let vertices;
let vertex_buffer;
let vertCode; 
let vertShader;
let vertStatus;
let vertErrors;
let fragCode; 
let fragShader;
let fragStatus;
let fragErrors;
let shaderProgram;
let coord;
let counter = 0;
let uTime = 0;
let uniformLocation_uTime;
let oldTime, newTime, frameTime;

function addLineNumbers(string)
{
	const lines = string.split('\n');

	for (let i = 0; i < lines.length; i++)
		lines[i] = (i + 1) + ': ' + lines[i];

	return lines.join('\n');
}

/* Step1: Prepare the canvas and get WebGL context */

// get the canvas element using the DOM
canvas = document.getElementById('mycanvas');

gl = canvas.getContext('webgl');


/* Step2: Define the geometry and store it in buffer objects */

vertices = [
	-1.0, -1.0, // left triangle
	1.0, -1.0,
	-1.0, 1.0,

	-1.0, 1.0, // right triangle
	1.0, -1.0,
	1.0, 1.0
];

// Create a new buffer object
vertex_buffer = gl.createBuffer();

// Bind an empty array buffer to it
gl.bindBuffer(gl.ARRAY_BUFFER, vertex_buffer);

// Pass the vertices data to the buffer
gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);

// Unbind the buffer
gl.bindBuffer(gl.ARRAY_BUFFER, null);


/* Step3: Create and compile Shader programs */

// Vertex shader source code
vertCode = `
precision highp float;
attribute vec2 coordinates;

void main(void)
{
	gl_Position = vec4(coordinates, 0.0, 1.0);
}
`;

//Create a vertex shader object
vertShader = gl.createShader(gl.VERTEX_SHADER);

//Attach vertex shader source code
gl.shaderSource(vertShader, vertCode);

//Compile the vertex shader
gl.compileShader(vertShader);

vertStatus = gl.getShaderParameter(vertShader, gl.COMPILE_STATUS);
vertErrors = gl.getShaderInfoLog(vertShader);
if (!vertStatus || vertErrors != '')
	console.error(vertErrors + '\n' + addLineNumbers(gl.getShaderSource(vertShader)));


//Fragment shader source code
fragCode = `
precision highp float;

#define PI 3.14159265358979323
#define INFINITY 1000000.0

uniform float uTime;

struct Ray { vec3 origin; vec3 direction; };

void solveQuadratic(float A, float B, float C, out float t0, out float t1) // required for scenes with quadric shapes (spheres, cylinders, etc.)
{
	float invA = 1.0 / A;
	B *= invA;
	C *= invA;
	float neg_halfB = -B * 0.5;
	float u2 = neg_halfB * neg_halfB - C;
	float u = u2 < 0.0 ? neg_halfB = 0.0 : sqrt(u2);
	t0 = neg_halfB - u;
	t1 = neg_halfB + u;
}

float SphereIntersect( float rad, vec3 pos, Ray ray )
{
	float t0, t1;
	vec3 L = ray.origin - pos;
	float a = dot( ray.direction, ray.direction );
	float b = 2.0 * dot( ray.direction, L );
	float c = dot( L, L ) - (rad * rad);
	solveQuadratic(a, b, c, t0, t1);
	return t0 > 0.0 ? t0 : t1 > 0.0 ? t1 : INFINITY;
}


float PosY_XZRectangleIntersect( vec3 pos, float radiusU, float radiusV, Ray r )
{
	vec3 normal = vec3(0,1,0);
	float dt = dot(-normal, r.direction);
	// use the following for one-sided rectangle
	//if (dt < 0.0) return INFINITY;
	float t = dot(-normal, pos - r.origin) / dt;
	if (t < 0.0) return INFINITY;

	vec3 hit = r.origin + r.direction * t;
	vec3 rectPosToHit = hit - pos;
	vec3 U = vec3(1,0,0);
	vec3 V = vec3(0,0,-1);
	return (abs(dot(rectPosToHit, U)) > radiusU || abs(dot(rectPosToHit, V)) > radiusV) ? INFINITY : t;
}

void main(void)
{
	vec3 worldUp = vec3(0, 1, 0);

	float sphereRadius = 3.0;
	vec3 spherePosition = vec3(0,sphereRadius,0);

	float angle = uTime * (PI / 180.0) * 30.0;
	angle = mod(angle, 2.0 * PI);
	vec3 cameraPosition = spherePosition + vec3(cos(angle) * 10.0, 10, sin(angle) * 10.0);
	//cameraPosition = vec3(0,6,15);
	
	vec3 camRight, camUp, camForward;
	vec3 lookTarget = spherePosition;
	vec3 lookDirection = normalize( lookTarget - cameraPosition );
	// camera looks down -Z axis, but the 'camForward' basis vector must point the opposite way, so that it points to +Z axis
	// in the end, we need all camera basis vectors pointing to +X, +Y, and +Z axes: 
	// camRight points along the +X axis, camUp points along the +Y axis, 'camForward' points along the +Z axis
	camForward = -lookDirection; // by flipping this around, we get back to the true 'camForward' basis vector, (even though the camera 'looks' in the opposite Z direction)
	
	camRight = cross(worldUp, camForward);
	camRight = normalize(camRight);
	camUp = cross(camForward, camRight);
	camUp = normalize(camUp);

	vec2 resolution = vec2(640, 480);
	float aspectRatio = resolution.x / resolution.y;
	float fov = 60.0;
	float vLen = tan(fov * 0.5 * (PI / 180.0));
	float uLen = vLen * aspectRatio;
	vec2 uv = gl_FragCoord.xy / resolution;
	vec2 pixelPos = uv * 2.0 - 1.0;
	vec3 rayDir = normalize( pixelPos.x * camRight * uLen + pixelPos.y * camUp * vLen + lookDirection );

	vec3 pixelColor = vec3(0);
	Ray ray = Ray(cameraPosition, rayDir);

	float t = INFINITY;
	float d;

	d = PosY_XZRectangleIntersect( vec3(0,0,0), 10.0, 10.0, ray );
	if (d < t)
	{
		t = d;
		pixelColor = vec3(1,0,0);
	}
	d = SphereIntersect( sphereRadius, spherePosition, ray );
	if (d < t)
	{
		t = d;
		pixelColor = vec3(1,1,1);
	}

	gl_FragColor = vec4(pixelColor, 1.0);
}
`;


// Create fragment shader object
fragShader = gl.createShader(gl.FRAGMENT_SHADER);

// Attach fragment shader source code
gl.shaderSource(fragShader, fragCode);

// Compile the fragment shader
gl.compileShader(fragShader);

fragStatus = gl.getShaderParameter(fragShader, gl.COMPILE_STATUS);
fragErrors = gl.getShaderInfoLog(fragShader);
if (!fragStatus || fragErrors != '')
	console.error(fragErrors + '\n' + addLineNumbers(gl.getShaderSource(fragShader)));


// Create a shader program object to store combined shader program
shaderProgram = gl.createProgram();

// Attach a vertex shader
gl.attachShader(shaderProgram, vertShader);

// Attach a fragment shader
gl.attachShader(shaderProgram, fragShader);

// Link both programs
gl.linkProgram(shaderProgram);

// Use the combined shader program object
gl.useProgram(shaderProgram);


/* Step 4: Associate the shader programs to buffer objects */

//Bind vertex buffer object
gl.bindBuffer(gl.ARRAY_BUFFER, vertex_buffer);

//Get the attribute location
coord = gl.getAttribLocation(shaderProgram, "coordinates");

//point an attribute to the currently bound VBO
gl.vertexAttribPointer(coord, 2, gl.FLOAT, false, 0, 0);

//Enable the attribute
gl.enableVertexAttribArray(coord);

// set up uniform locations
uniformLocation_uTime = gl.getUniformLocation(shaderProgram, 'uTime');



/* Step5: Drawing the required objects (2 triangles) */

function drawScreenQuad()
{
	// Clear the canvas
	gl.clearColor(0.5, 0.5, 0.5, 1.0);

	// Enable the depth test
	//gl.enable(gl.DEPTH_TEST);

	// Clear the color buffer bit
	gl.clear(gl.COLOR_BUFFER_BIT);

	// Set the view port
	gl.viewport(0, 0, canvas.width, canvas.height);

	// Draw the triangle
	//gl.drawArrays(gl.LINES, 0, 6);
	gl.drawArrays(gl.TRIANGLES, 0, 6);
}



function getTime()
{
	newTime = performance.now();
	frameTime = newTime - oldTime;
	frameTime *= 0.001;
	uTime += frameTime;

	oldTime = newTime;
}

function animate() 
{
	getTime();

	// refresh uniforms
	gl.uniform1f(uniformLocation_uTime, uTime);

	drawScreenQuad();
	//infoElement.innerHTML = "Samples: " + sampleCount;

	requestAnimationFrame(animate);
}

oldTime = performance.now();

// start up the progressive rendering!
animate();
