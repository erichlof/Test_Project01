let canvas;
let gl;
let stats;
let vertices;
let vertex_buffer;
let commonVertShader;
let vertStatus;
let vertErrors;
let rayTracingFragShader;
let rayTracingShaderProgram;
let screenCopyFragShader;
let screenCopyShaderProgram;
let screenOutputFragShader;
let screenOutputShaderProgram;
let fragStatus;
let fragErrors;
let coord;
let counter = 0;
let resolutionX, resolutionY;
let uTime = 0;
let uniformLocation_uRayTracedImageTexture;
let uniformLocation_uPreviousScreenImageTexture;
let uniformLocation_uAccumulationBufferTexture;
let uniformLocation_uBlueNoiseTexture;
let uniformLocation_uDiffuseTexture;
let uniformLocation_uResolution;
let uniformLocation_uRandomVec2;
let uniformLocation_uApertureSize;
let uniformLocation_uFocusDistance;
let uniformLocation_uTime;
let oldTime, newTime, frameTime;
let uFrameCounter = 1;
let uniformLocation_uFrameCounter;
let sampleCounter = 0;
let uOneOverSampleCounter = 1;
let uniformLocation_uOneOverSampleCounter;
let windowIsBeingResized = true;
let uSceneIsDynamic;
let uniformLocation_uSceneIsDynamic;
let mouseControl = true;
let isPaused = true;
let pointerlockChange;
let apertureSize = 0.0;
let increaseAperture = false;
let decreaseAperture = false;
let focusDistance = 10.0;
let increaseFocusDist = false;
let decreaseFocusDist = false;
let increaseFOV = false;
let decreaseFOV = false;
let aspectRatio;
let FOV;
let uULen, uVLen;
let uniformLocation_uULen, uniformLocation_uVLen;
let displayWidth, displayHeight;
let needResize;
let image1, image2;
let uBlueNoiseTexture;
let uDiffuseTexture;
let angle, scale;
let scene = new Object3D();
let worldCamera = new Object3D();
let sphere0 = new Object3D();
let sphere1 = new Object3D();
let box0 = new Object3D();
let box1 = new Object3D();
let shearMatrix = new Matrix4();
let uSphere0InvMatrix = new Matrix4();
let uniformLocation_uSphere0InvMatrix;
let uSphere1InvMatrix = new Matrix4();
let uniformLocation_uSphere1InvMatrix;
let uBox0InvMatrix = new Matrix4();
let uniformLocation_uBox0InvMatrix;
let uBox1InvMatrix = new Matrix4();
let uniformLocation_uBox1InvMatrix;
let controls, cameraControlsObject, cameraControlsYawObject, cameraControlsPitchObject;
let oldYawRotation;
let oldPitchRotation;
let button1Pressed = false, button2Pressed = false, button3Pressed = false,
	button4Pressed = false, button5Pressed = false, button6Pressed = false;
// the following variables will be used to calculate rotations and directions from the camera
let cameraDirectionVector = new Vector3(); //for moving where the camera is looking
let cameraRightVector = new Vector3(); //for strafing the camera right and left
let cameraUpVector = new Vector3(); //for moving camera up and down
let cameraWorldQuaternion = new Quaternion(); //for rotating scene objects to match camera's current rotation
let cameraRotationSpeed = 1;
let camFlightSpeed = 10;
let uCameraIsMoving = false;
let cameraRecentlyMoving = false;
let mobileInPortraitMode;
let uniformLocation_uCameraIsMoving;
let uniformLocation_uCameraMatrix;

let uMatrices = [];
const numOfMatrices = 2;
let uniformLocation_uMatrices;
let elArray;
let matricesElementsArray = new Float32Array(16 * numOfMatrices);
for (let i = 0; i < numOfMatrices; i++)
{
	uMatrices.push(new Matrix4());
}

let cameraInfoElement = document.getElementById('cameraInfo');
cameraInfoElement.style.cursor = "default";
cameraInfoElement.style.userSelect = "none";
cameraInfoElement.style.MozUserSelect = "none";


let FirstPersonCameraControls = function (camera)
{
	camera.rotation.set(0, 0, 0);

	let pitchObject = new Object3D();
	pitchObject.add(camera);

	let yawObject = new Object3D();
	yawObject.add(pitchObject);

	let movementX = 0;
	let movementY = 0;

	var onMouseMove = function (event)
	{
		if (isPaused)
			return;

		movementX = event.movementX || event.mozMovementX || 0;
		movementY = event.movementY || event.mozMovementY || 0;

		yawObject.rotation.y -= movementX * 0.0012 * cameraRotationSpeed;
		pitchObject.rotation.x -= movementY * 0.001 * cameraRotationSpeed;
		// clamp the camera's vertical movement (around the x-axis) to the scene's 'ceiling' and 'floor'
		pitchObject.rotation.x = Math.max(- PI_2, Math.min(PI_2, pitchObject.rotation.x));
	};

	document.addEventListener('mousemove', onMouseMove, false);

	this.getObject = function ()
	{
		return yawObject;
	};

	this.getYawObject = function ()
	{
		return yawObject;
	};

	this.getPitchObject = function ()
	{
		return pitchObject;
	};

	this.getDirection = function ()
	{
		const te = pitchObject.matrixWorld.elements;

		return function (v)
		{
			v.set(te[8], te[9], te[10]).negate();
			return v;
		};
	}();

	this.getUpVector = function ()
	{
		const te = pitchObject.matrixWorld.elements;

		return function (v)
		{
			v.set(te[4], te[5], te[6]);
			return v;
		};
	}();

	this.getRightVector = function ()
	{
		const te = pitchObject.matrixWorld.elements;

		return function (v)
		{
			v.set(te[0], te[1], te[2]);
			return v;
		};
	}();
};



// Scene Setup /////////////////////////////////////////////////////////////////////////////////////

uSceneIsDynamic = true;

scene.add(worldCamera);

controls = new FirstPersonCameraControls(worldCamera);

cameraControlsObject = controls.getObject();
cameraControlsYawObject = controls.getYawObject();
cameraControlsPitchObject = controls.getPitchObject();

scene.add(cameraControlsObject);

cameraControlsObject.position.set(0, 5, 10);
cameraControlsPitchObject.rotation.set(Math.PI * -0.2, 0, 0);

FOV = 60;



function addLineNumbers(string)
{
	const lines = string.split('\n');

	for (let i = 0; i < lines.length; i++)
		lines[i] = (i + 1) + ': ' + lines[i];

	return lines.join('\n');
}

const KEYCODE_NAMES = {
	65: 'a', 66: 'b', 67: 'c', 68: 'd', 69: 'e', 70: 'f', 71: 'g', 72: 'h', 73: 'i', 74: 'j', 75: 'k', 76: 'l', 77: 'm',
	78: 'n', 79: 'o', 80: 'p', 81: 'q', 82: 'r', 83: 's', 84: 't', 85: 'u', 86: 'v', 87: 'w', 88: 'x', 89: 'y', 90: 'z',
	37: 'left', 38: 'up', 39: 'right', 40: 'down', 32: 'space', 33: 'pageup', 34: 'pagedown', 9: 'tab',
	189: 'dash', 187: 'equals', 188: 'comma', 190: 'period', 27: 'escape', 13: 'enter'
}
let KeyboardState = {
	a: false, b: false, c: false, d: false, e: false, f: false, g: false, h: false, i: false, j: false, k: false, l: false, m: false,
	n: false, o: false, p: false, q: false, r: false, s: false, t: false, u: false, v: false, w: false, x: false, y: false, z: false,
	left: false, up: false, right: false, down: false, space: false, pageup: false, pagedown: false, tab: false,
	dash: false, equals: false, comma: false, period: false, escape: false, enter: false
}

function onKeyDown(event)
{
	event.preventDefault();

	KeyboardState[KEYCODE_NAMES[event.keyCode]] = true;
}

function onKeyUp(event)
{
	event.preventDefault();

	KeyboardState[KEYCODE_NAMES[event.keyCode]] = false;
}

function keyPressed(keyName)
{
	if (!mouseControl)
		return;

	return KeyboardState[keyName];
}

function onMouseWheel(event)
{
	if (isPaused)
		return;

	// use the following instead, because event.preventDefault() gives errors in console
	event.stopPropagation();

	if (event.deltaY > 0)
	{
		increaseFOV = true;
	}
	else if (event.deltaY < 0)
	{
		decreaseFOV = true;
	}
}


window.addEventListener('resize', onWindowResize, false);


if ('ontouchstart' in window)
{
	mouseControl = false;

	/* mobileJoystickControls = new MobileJoystickControls({
		//showJoystick: true
	}); */
}

// if on mobile device, unpause the app because there is no ESC key and no mouse capture to do
if (!mouseControl)
	isPaused = false;

if (mouseControl)
{
	window.addEventListener('wheel', onMouseWheel, false);

	document.body.addEventListener("click", function ()
	{
		this.requestPointerLock = this.requestPointerLock || this.mozRequestPointerLock;
		this.requestPointerLock();
	}, false);

	window.addEventListener("click", function (event)
	{
		event.preventDefault();
	}, false);

	window.addEventListener("dblclick", function (event)
	{
		event.preventDefault();
	}, false);


	pointerlockChange = function (event)
	{
		if (document.pointerLockElement === document.body ||
			document.mozPointerLockElement === document.body || document.webkitPointerLockElement === document.body)
		{
			document.addEventListener('keydown', onKeyDown, false);
			document.addEventListener('keyup', onKeyUp, false);
			isPaused = false;
		}
		else
		{
			document.removeEventListener('keydown', onKeyDown, false);
			document.removeEventListener('keyup', onKeyUp, false);
			isPaused = true;
		}
	};

	// Hook pointer lock state change events
	document.addEventListener('pointerlockchange', pointerlockChange, false);
	document.addEventListener('mozpointerlockchange', pointerlockChange, false);
	document.addEventListener('webkitpointerlockchange', pointerlockChange, false);
}

/*
// Fullscreen API
document.addEventListener("click", function() {
	
	if ( !document.fullscreenElement && !document.mozFullScreenElement && !document.webkitFullscreenElement ) 
	{
		if (document.documentElement.requestFullscreen) 
		{
			document.documentElement.requestFullscreen();	
		} 
		else if (document.documentElement.mozRequestFullScreen) 
		{
			document.documentElement.mozRequestFullScreen();
		} else if (document.documentElement.webkitRequestFullscreen) 
		{
			document.documentElement.webkitRequestFullscreen();
		}
	}
});
*/



/* Step1: Prepare the canvas and get WebGL context */

// get the canvas element using the DOM
canvas = document.getElementById('mycanvas');

gl = canvas.getContext('webgl2');

gl.getExtension('EXT_color_buffer_float');

stats = new Stats();
stats.dom.style.position = 'absolute';
stats.dom.style.top = '0px';
stats.dom.style.cursor = "default";
stats.dom.style.webkitUserSelect = "none";
stats.dom.style.MozUserSelect = "none";
document.body.appendChild(stats.dom);


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



// 
// Create a blue noise texture.
uBlueNoiseTexture = gl.createTexture();

// Asynchronously load an image
image1 = new Image();
image1.src = 'textures/BlueNoise_RGBA256.png';
image1.addEventListener('load', function ()
{
	// use texture unit 1
	gl.activeTexture(gl.TEXTURE1);
	// Now that the image has loaded, copy it to the texture.
	gl.bindTexture(gl.TEXTURE_2D, uBlueNoiseTexture);
	gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
	gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image1);
	//gl.generateMipmap(gl.TEXTURE_2D);
	gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
});

// Create a regular diffuse color texture.
uDiffuseTexture = gl.createTexture();

// Asynchronously load an image
image2 = new Image();
image2.src = 'textures/f-texture.png';
image2.addEventListener('load', function ()
{
	// use texture unit 2
	gl.activeTexture(gl.TEXTURE2);
	// Now that the image has loaded, copy it to the texture.
	gl.bindTexture(gl.TEXTURE_2D, uDiffuseTexture);
	gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
	gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image2);
	//gl.generateMipmap(gl.TEXTURE_2D);
	gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
});



// create a target texture to render the raytraced image to
const uRayTracedImageTexture = gl.createTexture();
gl.bindTexture(gl.TEXTURE_2D, uRayTracedImageTexture);
gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, 256, 256, 0, gl.RGBA, gl.FLOAT, null);
// set the filtering so we don't need mips
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

// Create and bind the framebuffer
const rayTracedImage_FrameBuffer = gl.createFramebuffer();
gl.bindFramebuffer(gl.FRAMEBUFFER, rayTracedImage_FrameBuffer);
// attach the texture as the first color attachment
gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, uRayTracedImageTexture, 0);


// create a target texture to render the copied screen image to
const uPreviousScreenImageTexture = gl.createTexture();
gl.bindTexture(gl.TEXTURE_2D, uPreviousScreenImageTexture);
gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, 256, 256, 0, gl.RGBA, gl.FLOAT, null);
// set the filtering so we don't need mips
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

// Create and bind the framebuffer
const previousScreenImage_FrameBuffer = gl.createFramebuffer();
gl.bindFramebuffer(gl.FRAMEBUFFER, previousScreenImage_FrameBuffer);
// attach the texture as the first color attachment
gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, uPreviousScreenImageTexture, 0);








/* Step3: Create and compile Shader programs */

//Create a vertex shader object
commonVertShader = gl.createShader(gl.VERTEX_SHADER);

//Attach vertex shader source code
gl.shaderSource(commonVertShader, commonVertCode);

//Compile the vertex shader
gl.compileShader(commonVertShader);

vertStatus = gl.getShaderParameter(commonVertShader, gl.COMPILE_STATUS);
vertErrors = gl.getShaderInfoLog(commonVertShader);
if (!vertStatus || vertErrors != '')
	console.error(vertErrors + '\n' + addLineNumbers(gl.getShaderSource(commonVertShader)));



////////////////////////////////////////////////////////////////////
// Create screenCopy shader object
screenCopyFragShader = gl.createShader(gl.FRAGMENT_SHADER);

// Attach fragment shader source code
gl.shaderSource(screenCopyFragShader, screenCopyFragCode);

// Compile the fragment shader
gl.compileShader(screenCopyFragShader);

fragStatus = gl.getShaderParameter(screenCopyFragShader, gl.COMPILE_STATUS);
fragErrors = gl.getShaderInfoLog(screenCopyFragShader);
if (!fragStatus || fragErrors != '')
	console.error(fragErrors + '\n' + addLineNumbers(gl.getShaderSource(screenCopyFragShader)));

// Create a shader program object to store combined shader program
screenCopyShaderProgram = gl.createProgram();

// Attach a vertex shader
gl.attachShader(screenCopyShaderProgram, commonVertShader);

// Attach a fragment shader
gl.attachShader(screenCopyShaderProgram, screenCopyFragShader);

// Link both programs
gl.linkProgram(screenCopyShaderProgram);



////////////////////////////////////////////////////////////////////
// Create final screenOutput shader object
screenOutputFragShader = gl.createShader(gl.FRAGMENT_SHADER);

// Attach fragment shader source code
gl.shaderSource(screenOutputFragShader, screenOutputFragCode);

// Compile the fragment shader
gl.compileShader(screenOutputFragShader);

fragStatus = gl.getShaderParameter(screenOutputFragShader, gl.COMPILE_STATUS);
fragErrors = gl.getShaderInfoLog(screenOutputFragShader);
if (!fragStatus || fragErrors != '')
	console.error(fragErrors + '\n' + addLineNumbers(gl.getShaderSource(screenOutputFragShader)));

// Create a shader program object to store combined shader program
screenOutputShaderProgram = gl.createProgram();

// Attach a vertex shader
gl.attachShader(screenOutputShaderProgram, commonVertShader);

// Attach a fragment shader
gl.attachShader(screenOutputShaderProgram, screenOutputFragShader);

// Link both programs
gl.linkProgram(screenOutputShaderProgram);



///////////////////////////////////////////////////
// Create rayTracing fragment shader object
rayTracingFragShader = gl.createShader(gl.FRAGMENT_SHADER);

// Attach fragment shader source code
gl.shaderSource(rayTracingFragShader, rayTracingFragCode);

// Compile the fragment shader
gl.compileShader(rayTracingFragShader);

fragStatus = gl.getShaderParameter(rayTracingFragShader, gl.COMPILE_STATUS);
fragErrors = gl.getShaderInfoLog(rayTracingFragShader);
if (!fragStatus || fragErrors != '')
	console.error(fragErrors + '\n' + addLineNumbers(gl.getShaderSource(rayTracingFragShader)));

// Create a shader program object to store combined shader program
rayTracingShaderProgram = gl.createProgram();

// Attach a vertex shader
gl.attachShader(rayTracingShaderProgram, commonVertShader);

// Attach a fragment shader
gl.attachShader(rayTracingShaderProgram, rayTracingFragShader);

// Link both programs
gl.linkProgram(rayTracingShaderProgram);
///////////////////////////////////////////////////////////////////////



/* Step 4: Associate the shader programs to buffer objects */

//Bind vertex buffer object
gl.bindBuffer(gl.ARRAY_BUFFER, vertex_buffer);

//Get the attribute location
coord = gl.getAttribLocation(rayTracingShaderProgram, "coordinates");

//point an attribute to the currently bound VBO
gl.vertexAttribPointer(coord, 2, gl.FLOAT, false, 0, 0);

//Enable the attribute
gl.enableVertexAttribArray(coord);

//get access to uniform locations
uniformLocation_uRayTracedImageTexture = gl.getUniformLocation(screenCopyShaderProgram, 'uRayTracedImageTexture');

uniformLocation_uAccumulationBufferTexture = gl.getUniformLocation(screenOutputShaderProgram, 'uAccumulationBufferTexture');
uniformLocation_uOneOverSampleCounter = gl.getUniformLocation(screenOutputShaderProgram, 'uOneOverSampleCounter');

uniformLocation_uPreviousScreenImageTexture = gl.getUniformLocation(rayTracingShaderProgram, 'uPreviousScreenImageTexture');
uniformLocation_uBlueNoiseTexture = gl.getUniformLocation(rayTracingShaderProgram, 'uBlueNoiseTexture');
uniformLocation_uDiffuseTexture = gl.getUniformLocation(rayTracingShaderProgram, 'uDiffuseTexture');
uniformLocation_uMatrices = gl.getUniformLocation(rayTracingShaderProgram, 'uMatrices');
uniformLocation_uCameraMatrix = gl.getUniformLocation(rayTracingShaderProgram, 'uCameraMatrix');
uniformLocation_uSphere0InvMatrix = gl.getUniformLocation(rayTracingShaderProgram, 'uSphere0InvMatrix');
uniformLocation_uSphere1InvMatrix = gl.getUniformLocation(rayTracingShaderProgram, 'uSphere1InvMatrix');
uniformLocation_uBox0InvMatrix = gl.getUniformLocation(rayTracingShaderProgram, 'uBox0InvMatrix');
uniformLocation_uBox1InvMatrix = gl.getUniformLocation(rayTracingShaderProgram, 'uBox1InvMatrix');
uniformLocation_uCameraIsMoving = gl.getUniformLocation(rayTracingShaderProgram, 'uCameraIsMoving');
uniformLocation_uSceneIsDynamic = gl.getUniformLocation(rayTracingShaderProgram, 'uSceneIsDynamic');
uniformLocation_uTime = gl.getUniformLocation(rayTracingShaderProgram, 'uTime');
uniformLocation_uApertureSize = gl.getUniformLocation(rayTracingShaderProgram, 'uApertureSize');
uniformLocation_uFocusDistance = gl.getUniformLocation(rayTracingShaderProgram, 'uFocusDistance');
uniformLocation_uFrameCounter = gl.getUniformLocation(rayTracingShaderProgram, 'uFrameCounter');
uniformLocation_uResolution = gl.getUniformLocation(rayTracingShaderProgram, 'uResolution');
uniformLocation_uRandomVec2 = gl.getUniformLocation(rayTracingShaderProgram, 'uRandomVec2');
uniformLocation_uULen = gl.getUniformLocation(rayTracingShaderProgram, 'uULen');
uniformLocation_uVLen = gl.getUniformLocation(rayTracingShaderProgram, 'uVLen');


/* Step5: Drawing the required objects (2 triangles) */

function drawScreenQuad()
{
	// Clear the canvas
	gl.clearColor(0.5, 0.5, 0.5, 1.0);

	// Enable the depth test
	//gl.enable(gl.DEPTH_TEST);

	// Clear the color buffer bit
	gl.clear(gl.COLOR_BUFFER_BIT);

	// Set the viewport
	gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

	// Draw the 2 side-by-side triangles (making up a large quad that stretches across the entire viewport)
	//gl.drawArrays(gl.LINES, 0, 6);
	gl.drawArrays(gl.TRIANGLES, 0, 6);
}

function onWindowResize(event)
{
	windowIsBeingResized = true;
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

	// reset this every animation frame
	uCameraIsMoving = false;

	// the following gl.useProgram needs to be here in order to set the correct uniforms if windowIsBeingResized == true
	// STEP 1: Perform RayTracing then render to the target: uRayTracedImageTexture
	gl.useProgram(rayTracingShaderProgram);

	if (mobileInPortraitMode)
	{
		if (gl.canvas.clientWidth > gl.canvas.clientHeight)
		{
			windowIsBeingResized = true;
		}
	}
	if (!mobileInPortraitMode)
	{
		if (gl.canvas.clientWidth < gl.canvas.clientHeight)
		{
			windowIsBeingResized = true;
		}
	}


	if (windowIsBeingResized)
	{
		if (gl.canvas.clientWidth > gl.canvas.clientHeight)
		{
			mobileInPortraitMode = false;
		}
		else mobileInPortraitMode = true;

		uCameraIsMoving = true;

		gl.uniform1f(uniformLocation_uSceneIsDynamic, uSceneIsDynamic);

		displayWidth = gl.canvas.clientWidth;
		displayHeight = gl.canvas.clientHeight;

		if (devicePixelRatio > 1)
		{
			if (mobileInPortraitMode)
			{
				displayWidth *= 0.4;
				displayHeight *= 0.4;
			}
			else
			{
				displayWidth *= 0.8;
				displayHeight *= 0.8;
			}
		}
		else
		{
			displayWidth *= 0.85;
			displayHeight *= 0.85;
		}

		// Make the canvas the same size
		gl.canvas.width = displayWidth;
		gl.canvas.height = displayHeight;

		aspectRatio = displayWidth / displayHeight;
		uVLen = Math.tan(degToRad(FOV * 0.5));
		uULen = uVLen * aspectRatio;
		gl.uniform1f(uniformLocation_uULen, uULen);
		gl.uniform1f(uniformLocation_uVLen, uVLen);

		gl.uniform2f(uniformLocation_uResolution, displayWidth, displayHeight);

		// resize render target textures as well
		gl.bindTexture(gl.TEXTURE_2D, uRayTracedImageTexture);
		gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, displayWidth, displayHeight, 0, gl.RGBA, gl.FLOAT, null);

		gl.bindTexture(gl.TEXTURE_2D, uPreviousScreenImageTexture);
		gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, displayWidth, displayHeight, 0, gl.RGBA, gl.FLOAT, null);

		windowIsBeingResized = false;
	}



	// check user controls
	if (mouseControl)
	{
		// movement detected
		if (oldYawRotation != cameraControlsYawObject.rotation.y ||
			oldPitchRotation != cameraControlsPitchObject.rotation.x)
		{
			uCameraIsMoving = true;
		}

		// save state for next frame
		oldYawRotation = cameraControlsYawObject.rotation.y;
		oldPitchRotation = cameraControlsPitchObject.rotation.x;

	} // end if (mouseControl)

	// this gives us a vector in the direction that the camera is pointing,
	// which will be useful for moving the camera 'forward' and shooting projectiles in that direction
	controls.getDirection(cameraDirectionVector);
	cameraDirectionVector.normalize();
	controls.getUpVector(cameraUpVector);
	cameraUpVector.normalize();
	controls.getRightVector(cameraRightVector);
	cameraRightVector.normalize();
	// the following gives us a rotation quaternion (4D vector), which will be useful for 
	// rotating scene objects to match the camera's rotation
	worldCamera.getWorldQuaternion(cameraWorldQuaternion);


	// USER INPUT ////////////////////////////////////////////////////////////////////////
	if (!isPaused)
	{
		if ((keyPressed('w') || button3Pressed) && !(keyPressed('s') || button4Pressed))
		{
			cameraControlsObject.position.add(cameraDirectionVector.multiplyScalar(camFlightSpeed * frameTime));
			uCameraIsMoving = true;
		}
		if ((keyPressed('s') || button4Pressed) && !(keyPressed('w') || button3Pressed))
		{
			cameraControlsObject.position.sub(cameraDirectionVector.multiplyScalar(camFlightSpeed * frameTime));
			uCameraIsMoving = true;
		}
		if ((keyPressed('a') || button1Pressed) && !(keyPressed('d') || button2Pressed))
		{
			cameraControlsObject.position.sub(cameraRightVector.multiplyScalar(camFlightSpeed * frameTime));
			uCameraIsMoving = true;
		}
		if ((keyPressed('d') || button2Pressed) && !(keyPressed('a') || button1Pressed))
		{
			cameraControlsObject.position.add(cameraRightVector.multiplyScalar(camFlightSpeed * frameTime));
			uCameraIsMoving = true;
		}
		if (keyPressed('q') && !keyPressed('z'))
		{
			cameraControlsObject.position.add(cameraUpVector.multiplyScalar(camFlightSpeed * frameTime));
			uCameraIsMoving = true;
		}
		if (keyPressed('z') && !keyPressed('q'))
		{
			cameraControlsObject.position.sub(cameraUpVector.multiplyScalar(camFlightSpeed * frameTime));
			uCameraIsMoving = true;
		}

		if ((keyPressed('up') || button5Pressed) && !(keyPressed('down') || button6Pressed))
		{
			increaseFocusDist = true;
		}
		if ((keyPressed('down') || button6Pressed) && !(keyPressed('up') || button5Pressed))
		{
			decreaseFocusDist = true;
		}
		if (keyPressed('right') && !keyPressed('left'))
		{
			increaseAperture = true;
		}
		if (keyPressed('left') && !keyPressed('right'))
		{
			decreaseAperture = true;
		}
	} // end if (!isPaused)

	if (increaseFOV)
	{
		FOV++;
		if (FOV > 150)
			FOV = 150;

		uVLen = Math.tan(degToRad(FOV * 0.5));
		uULen = uVLen * aspectRatio;
		gl.uniform1f(uniformLocation_uULen, uULen);
		gl.uniform1f(uniformLocation_uVLen, uVLen);

		uCameraIsMoving = true;
		increaseFOV = false;
	}
	if (decreaseFOV)
	{
		FOV--;
		if (FOV < 1)
			FOV = 1;

		uVLen = Math.tan(degToRad(FOV * 0.5));
		uULen = uVLen * aspectRatio;
		gl.uniform1f(uniformLocation_uULen, uULen);
		gl.uniform1f(uniformLocation_uVLen, uVLen);

		uCameraIsMoving = true;
		decreaseFOV = false;
	}

	if (increaseFocusDist)
	{
		focusDistance += 0.1;
		
		uCameraIsMoving = true;
		increaseFocusDist = false;
	}
	if (decreaseFocusDist)
	{
		focusDistance -= 0.1;
		if (focusDistance < 1.0)
			focusDistance = 1.0;
		
		uCameraIsMoving = true;
		decreaseFocusDist = false;
	}

	if (increaseAperture)
	{
		apertureSize += 0.01;
		if (apertureSize > 10.0)
			apertureSize = 10.0;
		
		uCameraIsMoving = true;
		increaseAperture = false;
	}
	if (decreaseAperture)
	{
		apertureSize -= 0.01;
		if (apertureSize < 0.0)
			apertureSize = 0.0;
		
		uCameraIsMoving = true;
		decreaseAperture = false;
	}

	if (uFrameCounter > 1000)
		uFrameCounter = 1;

	if (!uCameraIsMoving)
	{
		if (uSceneIsDynamic)
			sampleCounter = 1.0; // reset for continuous updating of image
		else sampleCounter += 1.0; // for progressive refinement of image

		uFrameCounter += 1.0;

		cameraRecentlyMoving = false;
	}

	if (uCameraIsMoving)
	{
		sampleCounter = 1.0;
		uFrameCounter += 1.0;

		if (!cameraRecentlyMoving)
		{
			uFrameCounter = 1.0;
			cameraRecentlyMoving = true;
		}
	}


	// refresh uniforms
	gl.uniform2f(uniformLocation_uRandomVec2, Math.random(), Math.random());
	gl.uniform1f(uniformLocation_uCameraIsMoving, uCameraIsMoving);
	gl.uniform1f(uniformLocation_uTime, uTime);
	gl.uniform1f(uniformLocation_uFrameCounter, uFrameCounter);
	gl.uniform1f(uniformLocation_uFocusDistance, focusDistance);
	gl.uniform1f(uniformLocation_uApertureSize, apertureSize);

	elArray = [];
	for (let i = 0; i < numOfMatrices; i++)
	{
		for (let j = 0; j < 16; j++)
			elArray.push(uMatrices[i].elements[j]);
	}
	matricesElementsArray.set(elArray);
	gl.uniformMatrix4fv(uniformLocation_uMatrices, false, matricesElementsArray);

	//angle = degToRad(uTime * 10);
	//angle = angle % TWO_PI;

	// set camera Position
	//worldCamera.position.set(Math.cos(angle) * 10, Math.sin(angle) * 5 + 5.5, Math.sin(angle) * 10);

	// set camera Rotation
	//worldCamera.lookAt(scene.position);
	//worldCamera.rotateY(Math.PI);


	// send camera's Matrix to the GPU
	gl.uniformMatrix4fv(uniformLocation_uCameraMatrix, false, worldCamera.matrixWorld.elements);


	// SPHERE0
	// the following clears out object's matrix from last frame
	sphere0.updateMatrixWorld(true);
	// now build up a series of transformations on the object (position, rotation, scale, shear)
	//sphere0.position.set(0, Math.abs(Math.sin(uTime)) * 2 + 1, 6);
	sphere0.position.set(0, 2, 6);

	//sphere0.rotation.set(0, uTime % TWO_PI, 0);

	//scale = Math.abs(Math.sin(uTime)) * 2 + 0.1;
	//sphere0.scale.set(1, scale, 1);

	//shearMatrix.makeShear(0.0, 0.0, 0.0, 0.0, Math.sin(uTime), 0.0); // (xy, xz,  yx, yz,  zx, zy)
	//sphere0.matrixWorld.multiply(shearMatrix);

	// finally, send Sphere0's Matrix Inverse to the GPU
	uSphere0InvMatrix.copy(sphere0.matrixWorld).invert();
	gl.uniformMatrix4fv(uniformLocation_uSphere0InvMatrix, false, uSphere0InvMatrix.elements);


	// SPHERE1
	// the following clears out object's matrix from last frame
	sphere1.updateMatrixWorld(true);
	// now build up a series of transformations on the object (position, rotation, scale, shear)
	//sphere1.position.set(Math.cos(uTime * 0.2) * 5, 1, Math.sin(uTime * 0.2) * 5);
	sphere1.position.set(5, 2, 6);

	//sphere1.rotation.set(0, uTime % TWO_PI, 0);

	//scale = Math.abs(Math.sin(uTime)) * 2 + 0.1;
	//sphere1.scale.set(1, scale, 1);

	//shearMatrix.makeShear(0.0, 0.0, 0.0, 0.0, Math.sin(uTime), 0.0); // (xy, xz,  yx, yz,  zx, zy)
	//sphere1.matrixWorld.multiply(shearMatrix);

	// finally, send Sphere1's Matrix Inverse to the GPU
	uSphere1InvMatrix.copy(sphere1.matrixWorld).invert();
	gl.uniformMatrix4fv(uniformLocation_uSphere1InvMatrix, false, uSphere1InvMatrix.elements);


	// BOX0
	// the following clears out object's matrix from last frame
	box0.updateMatrixWorld(true);
	// now build up a series of transformations on the object (position, rotation, scale, shear)
	//box0.position.set(Math.sin(uTime) * 2, 1, 0);
	box0.position.set(-3, 2, 0);

	//box0.rotation.set(0, Math.PI * 0.25, Math.PI * 0.25);
	box0.rotation.set(0, Math.PI * 0.25, 0);

	//scale = Math.abs(Math.sin(uTime)) * 2 + 0.1;
	//scale = 1;
	//box0.scale.set(scale, scale, scale);

	//shearMatrix.makeShear(-0.1, 0.0, -0.1, 0.0, 0.0, 0.0); // (xy, xz, yx, yz, zx, zy)
	//box0.matrixWorld.multiply(shearMatrix);

	// finally, send box0's Matrix Inverse to the GPU
	uBox0InvMatrix.copy(box0.matrixWorld).invert();
	gl.uniformMatrix4fv(uniformLocation_uBox0InvMatrix, false, uBox0InvMatrix.elements);


	// BOX1
	// the following clears out object's matrix from last frame
	box1.updateMatrixWorld(true);
	// now build up a series of transformations on the object (position, rotation, scale, shear)
	//box1.position.set(Math.sin(uTime) * 2, 1, 0);
	box1.position.set(-3, 2, 0);

	box1.rotation.set(Math.PI * 0.25, 0, 0);

	//scale = Math.abs(Math.sin(uTime)) * 2 + 0.1;
	//scale = 0.25;
	//box1.scale.set(1, 1, scale);

	//shearMatrix.makeShear(-0.1, 0.0, -0.1, 0.0, 0.0, 0.0); // (xy, xz, yx, yz, zx, zy)
	//box1.matrixWorld.multiply(shearMatrix);

	// finally, send box1's Matrix Inverse to the GPU
	uBox1InvMatrix.copy(box1.matrixWorld).invert();
	gl.uniformMatrix4fv(uniformLocation_uBox1InvMatrix, false, uBox1InvMatrix.elements);



	gl.activeTexture(gl.TEXTURE0);
	gl.bindTexture(gl.TEXTURE_2D, uPreviousScreenImageTexture);
	// Tell the shader to use texture unit 0 for uPreviousScreenImageTexture
	gl.uniform1i(uniformLocation_uPreviousScreenImageTexture, 0);

	gl.activeTexture(gl.TEXTURE1);
	gl.bindTexture(gl.TEXTURE_2D, uBlueNoiseTexture);
	// Tell the shader to use texture unit 1 for uBlueNoiseTexture
	gl.uniform1i(uniformLocation_uBlueNoiseTexture, 1);


	gl.activeTexture(gl.TEXTURE2);
	gl.bindTexture(gl.TEXTURE_2D, uDiffuseTexture);
	// Tell the shader to use texture unit 2 for uDiffuseTexture
	gl.uniform1i(uniformLocation_uDiffuseTexture, 2);

	// render to the targetTexture
	gl.bindFramebuffer(gl.FRAMEBUFFER, rayTracedImage_FrameBuffer);
	drawScreenQuad();


	// STEP 2: Screen Copy
	gl.useProgram(screenCopyShaderProgram);

	gl.activeTexture(gl.TEXTURE0);
	// copy the scene with the texture we just rendered to
	gl.bindTexture(gl.TEXTURE_2D, uRayTracedImageTexture);
	// Tell the shader to use texture unit 0 for uRayTracedImageTexture
	gl.uniform1i(uniformLocation_uRayTracedImageTexture, 0);

	// render to the target texture
	gl.bindFramebuffer(gl.FRAMEBUFFER, previousScreenImage_FrameBuffer);
	drawScreenQuad();


	// STEP 3: Screen Output
	gl.useProgram(screenOutputShaderProgram);

	gl.activeTexture(gl.TEXTURE0);
	// copy the scene with the texture we just rendered to
	gl.bindTexture(gl.TEXTURE_2D, uRayTracedImageTexture);
	// Tell the shader to use texture unit 0 for uRayTracedImageTexture
	gl.uniform1i(uniformLocation_uAccumulationBufferTexture, 0);


	uOneOverSampleCounter = 1.0 / sampleCounter;
	gl.uniform1f(uniformLocation_uOneOverSampleCounter, uOneOverSampleCounter);

	// render to the canvas
	gl.bindFramebuffer(gl.FRAMEBUFFER, null);
	drawScreenQuad();


	//cameraInfoElement.innerHTML = "FOV: " + FOV + "<br>" + "Samples: " + sampleCounter;
	//cameraInfoElement.innerHTML = "canvasWidth: " + document.body.clientWidth + " canvasHeight: " + document.body.clientHeight + "<br>" + "portraitMode: " + mobileInPortraitMode;
	cameraInfoElement.innerHTML = "displayWidth: " + displayWidth + " displayHeight: " + displayHeight + "<br>" + "portraitMode: " + mobileInPortraitMode;

	stats.update();

	requestAnimationFrame(animate);
}

oldTime = performance.now();

// start up the progressive rendering!
animate();