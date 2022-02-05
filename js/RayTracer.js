let canvas;
let gl;
let stats;
let vertices;
let vertex_buffer;
let commonVertCode;
let commonVertShader;
let vertStatus;
let vertErrors;
let rayTracingFragCode;
let rayTracingFragShader;
let rayTracingShaderProgram;
let screenCopyFragCode;
let screenCopyFragShader;
let screenCopyShaderProgram;
let screenOutputFragCode;
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
let uniformLocation_uDiffuseTexture;
let uniformLocation_uResolution;
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
let increaseFOV = false;
let decreaseFOV = false;
let aspectRatio;
let FOV;
let uULen, uVLen;
let uniformLocation_uULen, uniformLocation_uVLen;
let displayWidth, displayHeight;
let needResize;
let image;
let uDiffuseTexture;
let angle, scale;
let scene = new Object3D();
let worldCamera = new Object3D();
let sphere0 = new Object3D();
let box0 = new Object3D();
let shearMatrix = new Matrix4();
let uSphere0InvMatrix = new Matrix4();
let uniformLocation_uSphere0InvMatrix;
let uBox0InvMatrix = new Matrix4();
let uniformLocation_uBox0InvMatrix;
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


let FirstPersonCameraControls = function(camera)
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




// Create a regular diffuse color texture.
uDiffuseTexture = gl.createTexture();

// Asynchronously load an image
image = new Image();
image.src = "textures/f-texture.png";
image.addEventListener('load', function ()
{
	// use texture unit 1
	gl.activeTexture(gl.TEXTURE1);
	// Now that the image has loaded, copy it to the texture.
	gl.bindTexture(gl.TEXTURE_2D, uDiffuseTexture);
	gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
	gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);
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

// common Vertex shader source code
commonVertCode =
`#version 300 es

precision highp float;
precision highp int;
precision highp sampler2D;

in vec2 coordinates;

void main(void)
{
	gl_Position = vec4(coordinates, 0.0, 1.0);
}
`;

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

//screenCopyFragShader source code
screenCopyFragCode =
`#version 300 es

precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D uRayTracedImageTexture;

out vec4 pixelColor; // final pixel color output

void main()
{	
	pixelColor = texelFetch(uRayTracedImageTexture, ivec2(gl_FragCoord.xy), 0);
}
`;

//screenOutputFragShader source code
screenOutputFragCode =
`#version 300 es

precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D uAccumulationBufferTexture;
uniform float uOneOverSampleCounter;

out vec4 pixelColor; // final pixel color output

vec3 ReinhardToneMapping(vec3 color)
{
	///color *= uToneMappingExposure;
	return clamp(color / (vec3(1) + color), 0.0, 1.0);
}

void main()
{	
	pixelColor = texelFetch(uAccumulationBufferTexture, ivec2(gl_FragCoord.xy), 0);

	pixelColor *= uOneOverSampleCounter;

	// apply tone mapping (brings pixel into 0.0-1.0 rgb color range)
	pixelColor.rgb = ReinhardToneMapping(pixelColor.rgb);

	// lastly, apply gamma correction (gives more intensity/brightness range where it's needed)
	pixelColor = clamp(vec4( pow(pixelColor.rgb, vec3(0.4545)), 1.0 ), 0.0, 1.0);
}
`;


//rayTracingFragShader source code
rayTracingFragCode =
`#version 300 es

precision highp float;
precision highp int;
precision highp sampler2D;

#define PI 3.14159265358979323
#define TWO_PI 6.28318530717958648
#define INFINITY 1000000.0
#define SPOT_LIGHT -3
#define POINT_LIGHT -2
#define DIRECTIONAL_LIGHT -1
#define CHECKER 1
#define DIFFUSE 2
#define METAL 3
#define TRANSPARENT 4
#define CLEARCOAT_DIFFUSE 5

#define N_DIRECTIONAL_LIGHTS 1
#define N_POINT_LIGHTS 1
#define N_SPOT_LIGHTS 1

uniform sampler2D uPreviousScreenImageTexture;
uniform sampler2D uDiffuseTexture;
uniform mat4 uMatrices[2];
uniform mat4 uCameraMatrix;
uniform mat4 uSphere0InvMatrix;
uniform mat4 uBox0InvMatrix;
uniform vec2 uResolution;
uniform float uTime;
uniform float uFrameCounter;
uniform float uULen;
uniform float uVLen;
uniform bool uSceneIsDynamic;
uniform bool uCameraIsMoving;

out vec4 pixelColor; // final pixel color output

vec3 rayOrigin; 
vec3 rayDirection;
uvec2 seed;

struct DirectionalLight{ vec3 directionTowardsLight; vec3 color; float power; int type; };
struct PointLight{ vec3 position; vec3 color; float power; int type; };
struct SpotLight{ vec3 position; vec3 spotLightAimingDirection; float cutoffAngle; vec3 color; float power; int type; };

DirectionalLight directionalLights[N_DIRECTIONAL_LIGHTS];
PointLight pointLights[N_POINT_LIGHTS];
SpotLight spotLights[N_SPOT_LIGHTS];


vec3 ReinhardToneMapping(vec3 color)
{
	///color *= uToneMappingExposure;
	return clamp(color / (vec3(1) + color), 0.0, 1.0);
}

float calcFresnelReflectance(vec3 rayDirection, vec3 n, float etai, float etat, out float ratioIoR)
{
	float temp = etai;
	float cosi = clamp(dot(rayDirection, n), -1.0, 1.0);
	if (cosi > 0.0)
	{
		etai = etat;
		etat = temp;
	}
	
	ratioIoR = etai / etat;
	float sint = ratioIoR * sqrt(1.0 - (cosi * cosi));
	if (sint >= 1.0) 
		return 1.0; // total internal reflection
	float cost = sqrt(1.0 - (sint * sint));
	cosi = abs(cosi);
	float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
	float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
	return clamp( ((Rs * Rs) + (Rp * Rp)) * 0.5, 0.0, 1.0 );
}

float rng()
{
	seed += uvec2(1);
    	uvec2 q = 1103515245U * ( (seed >> 1U) ^ (seed.yx) );
    	uint  n = 1103515245U * ( (q.x) ^ (q.y >> 3U) );
	return float(n) * (1.0 / float(0xffffffffU));
}

/* vec3 randomDirectionInPhongSpecular(vec3 reflectionDir, float exponent)
{
	float phi = rng() * 2.0 * PI;
	float cosTheta = pow(rng(), 1.0 / (exponent + 1.0));
	float sinTheta = sqrt(max(0.0, 1.0 - cosTheta * cosTheta));
	
	vec3 U = normalize( cross( abs(reflectionDir.y) < 0.9 ? vec3(0, 1, 0) : vec3(1, 0, 0), reflectionDir ) );
	vec3 V = cross(reflectionDir, U);
	return normalize(U * cos(phi) * sinTheta + V * sin(phi) * sinTheta + reflectionDir * cosTheta);
} */

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

float SphereIntersect( float rad, vec3 pos, vec3 rayOrigin, vec3 rayDirection )
{
	float t0, t1;
	vec3 L = rayOrigin - pos;
	float a = dot( rayDirection, rayDirection );
	float b = 2.0 * dot( rayDirection, L );
	float c = dot( L, L ) - (rad * rad);
	solveQuadratic(a, b, c, t0, t1);
	return t0 > 0.0 ? t0 : t1 > 0.0 ? t1 : INFINITY;
}

//----------------------------------------------------------------------------------------------------------------------
float EllipsoidParamIntersect( float yMinPercent, float yMaxPercent, float phiMaxRadians, vec3 ro, vec3 rd, out vec3 n )
//----------------------------------------------------------------------------------------------------------------------
{
	vec3 pHit;
	float result = INFINITY; 
	float t0, t1, phi;
	// implicit equation of a unit (radius of 1) sphere:
	// x^2 + y^2 + z^2 - 1 = 0
	float a = dot(rd, rd);
	float b = 2.0 * dot(rd, ro);
	float c = dot(ro, ro) - 1.0;
	solveQuadratic(a, b, c, t0, t1);
	if (t1 > 0.0)
	{
		pHit = ro + rd * t1;
		phi = mod(atan(pHit.z, pHit.x), TWO_PI);
		if (pHit.y <= yMaxPercent && pHit.y >= yMinPercent && phi <= phiMaxRadians)
		{
			result = t1;
			n = vec3(2.0 * pHit.x, 2.0 * pHit.y, 2.0 * pHit.z);
		}	
	}
	if (t0 > 0.0)
	{
		pHit = ro + rd * t0;
		phi = mod(atan(pHit.z, pHit.x), TWO_PI);
		if (pHit.y <= yMaxPercent && pHit.y >= yMinPercent && phi <= phiMaxRadians)
		{
			result = t0;
			n = vec3(2.0 * pHit.x, 2.0 * pHit.y, 2.0 * pHit.z);
		}	
	}
	
	n = dot(rd, n) < 0.0 ? n : -n; // flip normal if it is facing away from us
	return result;
}

float PosY_XZRectangleIntersect( vec3 pos, float radiusU, float radiusV, vec3 rayOrigin, vec3 rayDirection )
{
	vec3 normal = vec3(0,1,0);
	float dt = dot(-normal, rayDirection);
	// use the following for one-sided rectangle
	//if (dt < 0.0) return INFINITY;
	float t = dot(-normal, pos - rayOrigin) / dt;
	if (t < 0.0) return INFINITY;

	vec3 hit = rayOrigin + rayDirection * t;
	vec3 rectPosToHit = hit - pos;
	vec3 U = vec3(1,0,0);
	vec3 V = vec3(0,0,-1);
	return (abs(dot(rectPosToHit, U)) > radiusU || abs(dot(rectPosToHit, V)) > radiusV) ? INFINITY : t;
}

float UnitBoxIntersect( vec3 ro, vec3 rd, out vec3 n )
{
	vec3 invDir = 1.0 / rd;
	vec3 near = (vec3(-1) - ro) * invDir; // unit radius box: vec3(-1,-1,-1) min corner
	vec3 far  = (vec3( 1) - ro) * invDir; // unit radius box: vec3(+1,+1,+1) max corner
	vec3 tmin = min(near, far);
	vec3 tmax = max(near, far);
	float t0 = max( max(tmin.x, tmin.y), tmin.z);
	float t1 = min( min(tmax.x, tmax.y), tmax.z);

        if (t0 < t1)
        {
                if (t0 > 0.0)
                {
                        n = -sign(rd) * step(tmin.yzx, tmin) * step(tmin.zxy, tmin);
                        return t0;
                }
                if (t1 > 0.0)
                {
                        n = -sign(rd) * step(tmax, tmax.yzx) * step(tmax, tmax.zxy);
                        return t1;
                }
        }

        return INFINITY;
}

float BoundingBoxIntersect( vec3 minCorner, vec3 maxCorner, vec3 rayOrigin, vec3 invDir )
{
	vec3 near = (minCorner - rayOrigin) * invDir;
	vec3 far  = (maxCorner - rayOrigin) * invDir;

	vec3 tmin = min(near, far);
	vec3 tmax = max(near, far);

	float t0 = max( max(tmin.x, tmin.y), tmin.z);
	float t1 = min( min(tmax.x, tmax.y), tmax.z);

	return max(t0, 0.0) > t1 ? INFINITY : t0;
}



float SceneIntersect(vec3 rayOrigin, vec3 rayDirection, bool shadowRayAimedAtPositionalLight, 
			out vec3 hitNormal, out vec3 hitColor, out float hitShininess, out int hitType)
{
	vec3 hitPos, n;
	vec3 rObjOrigin, rObjDirection;
	float t = INFINITY;
	float d;
	bool isRayExiting;

	d = PosY_XZRectangleIntersect( vec3(0,0,0), 10.0, 10.0, rayOrigin, rayDirection );
	if (d < t)
	{
		t = d;
		hitNormal = vec3(0, 1, 0);
		hitColor = vec3(1.0, 0.0, 0.0);
		hitShininess = 1000.0;
		hitType = CHECKER;
	}

	// transform ray into Ellipsoid Param's object space
	rObjOrigin = vec3( uSphere0InvMatrix * vec4(rayOrigin, 1.0) );
	rObjDirection = vec3( uSphere0InvMatrix * vec4(rayDirection, 0.0) );
	
	float angleAmount = (sin(uTime) * 0.5 + 0.5);
	//d = EllipsoidParamIntersect(-1.0, 1.0, TWO_PI, rObjOrigin, rObjDirection, n);
	d = EllipsoidParamIntersect(-1.0, angleAmount, PI, rObjOrigin, rObjDirection, n);
	if (d < t)
	{
		t = d;
		hitNormal = normalize(transpose(mat3(uSphere0InvMatrix)) * n);
		hitColor = vec3(0.0, 0.3, 1.0);
		hitShininess = 1000.0;
		hitType = CLEARCOAT_DIFFUSE;
	}

	// transform ray into Unit Box's object space
	rObjOrigin = vec3( uBox0InvMatrix * vec4(rayOrigin, 1.0) );
	rObjDirection = vec3( uBox0InvMatrix * vec4(rayDirection, 0.0) );
	
	d = UnitBoxIntersect(rObjOrigin, rObjDirection, n);
	if (d < t)
	{
		t = d;
		hitNormal = normalize(transpose(mat3(uBox0InvMatrix)) * n);
		hitColor = vec3(0.0, 1.0, 0.0);
		hitShininess = 1000.0;
		hitType = CLEARCOAT_DIFFUSE;
	}

	if (!shadowRayAimedAtPositionalLight)
		return t;

	d = BoundingBoxIntersect(pointLights[0].position - vec3(0.1), pointLights[0].position + vec3(0.1), rayOrigin, 1.0/rayDirection);
	if (d < t)
	{
		t = d;
		hitColor = pointLights[0].color * pointLights[0].power;
		hitType = POINT_LIGHT;
	}

	// float angleBetweenRayAndLight = acos(dot(rayDirection, -spotLights[0].spotLightAimingDirection));
	// if (angleBetweenRayAndLight > spotLights[0].cutoffAngle)
	// 	return t;

	// d = BoundingBoxIntersect(spotLights[0].position - vec3(0.1, 0.1, 0.25), spotLights[0].position + vec3(0.1, 0.1, 0.25), rayOrigin, 1.0/rayDirection);
	// if (d < t)
	// {
	// 	t = d;
	// 	hitColor = spotLights[0].color * spotLights[0].power;
	// 	hitType = SPOT_LIGHT;
	// }

	return t;
} // end float SceneIntersect()



vec3 getSkyColor(vec3 rayDirection)
{
	return mix( vec3(0), vec3(0.1, 0.8, 1.0), clamp((rayDirection.y + 1.0), 0.0, 1.0) );
}


vec3 CalculateRadiance()
{
	vec3 hitNormal, hitColor;
	float hitShininess;
	int hitType;

	vec3 checkColor0 = vec3(0.1);
	vec3 checkColor1 = vec3(1.0);
	vec3 colorMask = vec3(1);
	vec3 firstColorMask = vec3(1);
	vec3 finalColor = vec3(0);
	vec3 firstColor = vec3(0);
	vec3 ambientColor = vec3(0);
	vec3 diffuseColor = vec3(0);
	vec3 specularColor = vec3(0);
	vec3 n, nl, hitPoint;
	vec3 dirToLight;
	vec3 halfDirection;
	vec3 tdir;
	vec3 firstRayOrigin, firstRayDirection;

	float diffuseCosWeight = 0.0;
	float specularCosWeight = 0.0;
	float firstSpecularCosWeight = 0.0;
	float checkScale = 1.0;
	float t = INFINITY;
	float nc, nt, ratioIoR, Re, Tr;
	float P, RP, TP;

	bool bounceIsSpecular = true;
	bool isShadowRay = false;
	bool shadowRayAimedAtPositionalLight = false;
	bool firstHitWasTransparent = false;
	bool firstHitWasMetal = false;
	bool reflectionTime = false;

	for (int bounces = 0; bounces < 6; bounces++)
	{
		t = SceneIntersect(rayOrigin, rayDirection, shadowRayAimedAtPositionalLight, hitNormal, hitColor, hitShininess, hitType);
		
		if (t == INFINITY && !shadowRayAimedAtPositionalLight)
		{

			if (bounceIsSpecular)
			{
				finalColor = colorMask * getSkyColor(rayDirection);

				specularColor = colorMask * specularCosWeight * directionalLights[0].color * directionalLights[0].power;

				finalColor += specularColor;
			}
			if (isShadowRay) // shadow ray aimed towards light is not blocked by any scene object
			{
				// set object's ambient lighting (this is the darkest it can be, even in shadow)
				ambientColor = colorMask * 0.2;

				diffuseColor = colorMask * directionalLights[0].color * directionalLights[0].power;
				diffuseColor = mix(ambientColor, diffuseColor, diffuseCosWeight);

				specularColor = specularCosWeight * directionalLights[0].color * directionalLights[0].power;

				if (firstHitWasMetal)
					specularColor *= firstColorMask;
				finalColor = diffuseColor + specularColor; // diffuse color already includes ambient color in this method
			}

			if (firstHitWasTransparent)
			{
				if (!reflectionTime)
				{
					firstColor = finalColor;
					colorMask = firstColorMask;
					specularCosWeight = firstSpecularCosWeight;
					reflectionTime = true;
					bounceIsSpecular = true;
					isShadowRay = false;
					shadowRayAimedAtPositionalLight = false;
					rayDirection = firstRayDirection;
					rayOrigin = firstRayOrigin;
					continue;
				}
				else
				{
					finalColor += firstColor;
					break;
				}
			}

			break;
		} // end if (t == INFINITY)


		if (hitType == POINT_LIGHT || hitType == SPOT_LIGHT)
		{

			// set object's ambient lighting (this is the darkest it can be, even in shadow)
			ambientColor = colorMask * 0.2;

			diffuseColor = colorMask * hitColor;
			diffuseColor = mix(ambientColor, diffuseColor, diffuseCosWeight);

			specularColor = specularCosWeight * hitColor;

			if (firstHitWasMetal)
				specularColor *= firstColorMask;
			finalColor = diffuseColor + specularColor; // diffuse color already includes ambient color in this method
			
			if (firstHitWasTransparent)
			{
				if (!reflectionTime)
				{
					firstColor = finalColor;
					colorMask = firstColorMask;
					specularCosWeight = firstSpecularCosWeight;
					reflectionTime = true;
					bounceIsSpecular = true;
					isShadowRay = false;
					shadowRayAimedAtPositionalLight = false;
					rayDirection = firstRayDirection;
					rayOrigin = firstRayOrigin;
					continue;
				}
				else
				{
					finalColor += firstColor;
					break;
				}
			}

			break;
		} // end if (hitType == POINT_LIGHT || hitType == SPOT_LIGHT)

		
		// if still is shadow ray, the ray aimed towards the light is blocked by a scene object (it hit something), 
		// and therefore the surface remains in shadow 
		if (isShadowRay)
		{
			// ambientColor is for darker parts of the scene that are unlit by direct light sources
			// set object's ambient lighting (this is the darkest it can be, even in shadow)
			ambientColor = colorMask * 0.2;

			finalColor = ambientColor;

			if (firstHitWasTransparent)
			{
				if (!reflectionTime)
				{
					firstColor = finalColor;
					colorMask = firstColorMask;
					specularCosWeight = firstSpecularCosWeight;
					reflectionTime = true;
					bounceIsSpecular = true;
					isShadowRay = false;
					shadowRayAimedAtPositionalLight = false;
					rayDirection = firstRayDirection;
					rayOrigin = firstRayOrigin;
					continue;
				}
				else
				{
					finalColor += firstColor;
					break;
				}
			}

			break;
		}
			

		// useful data
		n = normalize(hitNormal);
		nl = dot(n, rayDirection) < 0.0 ? n : -n;
		hitPoint = rayOrigin + rayDirection * t;
		
		dirToLight = normalize(directionalLights[0].directionTowardsLight);
		//dirToLight= normalize(pointLights[0].position - hitPoint);
		//dirToLight= normalize(spotLights[0].position - hitPoint);


		if (hitType == DIFFUSE)
		{
			bounceIsSpecular = false;
			colorMask *= hitColor;

			// pre-calculate Lambert(cosine-weighted) diffuse lighting amount
			diffuseCosWeight = max(0.0, dot(dirToLight, nl));

			if (hitShininess > 0.0)
			{
				// pre-calculate Blinn-Phong(powered cosine) specular lighting amount
				halfDirection = normalize(-rayDirection + dirToLight);
				specularCosWeight = max(0.0, dot(halfDirection, nl));
				specularCosWeight = pow(specularCosWeight, hitShininess);
			}
			
			// create shadow ray
			rayDirection = dirToLight;
			rayOrigin = hitPoint + nl * 0.01;

			isShadowRay = true;
			//shadowRayAimedAtPositionalLight = true;
			continue;
		} // end if (hitType == DIFFUSE)


		if (hitType == METAL)
		{
			colorMask *= hitColor;

			if (bounces == 0)
			{
				firstHitWasMetal = true;
				firstColorMask = colorMask;
			}

			if (hitShininess > 0.0)
			{
				// pre-calculate Blinn-Phong(powered cosine) specular lighting amount
				halfDirection = normalize(-rayDirection + dirToLight);
				specularCosWeight = max(0.0, dot(halfDirection, nl));
				specularCosWeight = pow(specularCosWeight, hitShininess);
			}
			
			// create reflection ray
			rayDirection = reflect(rayDirection, nl);
			rayOrigin = hitPoint + nl * 0.01;

			continue;
		} // end if (hitType == METAL)


		if (hitType == TRANSPARENT)
		{
			nc = 1.0; // IOR of Air
			nt = 1.5; // IOR of common Glass
			Re = calcFresnelReflectance(rayDirection, n, nc, nt, ratioIoR);
			Tr = 1.0 - Re;
			P  = 0.25 + (0.5 * Re);
			RP = Re / P;
			TP = Tr / (1.0 - P);

			if (bounces == 0)
			{
				if (hitShininess > 0.0)
				{
					// pre-calculate Blinn-Phong(powered cosine) specular lighting amount
					halfDirection = normalize(-rayDirection + dirToLight);
					firstSpecularCosWeight = max(0.0, dot(halfDirection, nl));
					firstSpecularCosWeight = pow(firstSpecularCosWeight, hitShininess);
				}

				firstHitWasTransparent = true;
				firstColorMask *= Re;
				firstRayDirection = reflect(rayDirection, nl); // reflect ray from surface
				firstRayOrigin = hitPoint + nl * 0.01;
			}

			// transmit ray through surface
			colorMask *= hitColor;
			colorMask *= Tr;
			
			tdir = refract(rayDirection, nl, ratioIoR);
			rayDirection = tdir;
			rayOrigin = hitPoint - nl * 0.01;

			continue;
			
		} // end if (hitType == TRANSPARENT)

		if (hitType == CLEARCOAT_DIFFUSE || hitType == CHECKER)
		{
			nc = 1.0; // IOR of Air
			nt = 1.5; // IOR of clearCoat
			Re = calcFresnelReflectance(rayDirection, nl, nc, nt, ratioIoR);
			Tr = 1.0 - Re;
			P  = 0.25 + (0.5 * Re);
			RP = Re / P;
			TP = Tr / (1.0 - P);

			if (bounces == 0)
			{
				if (hitShininess > 0.0)
				{
					// pre-calculate Blinn-Phong(powered cosine) specular lighting amount
					halfDirection = normalize(-rayDirection + dirToLight);
					firstSpecularCosWeight = max(0.0, dot(halfDirection, nl));
					firstSpecularCosWeight = pow(firstSpecularCosWeight, hitShininess);
				}

				firstHitWasTransparent = true;
				firstColorMask *= Re;
				firstRayDirection = reflect(rayDirection, nl); // reflect ray from surface
				firstRayOrigin = hitPoint + nl * 0.01;
			}

			if (hitType == CHECKER)
			{
				//if (sin(hitPoint.x * checkScale) > -0.98 && sin(hitPoint.z * checkScale) > -0.98)
				/* if (mod(floor(hitPoint.x * checkScale), 2.0) + mod(floor(hitPoint.z * checkScale), 2.0) == 1.0)
					hitColor = checkColor0;
				else hitColor = checkColor1; */

				hitColor = pow(texture(uDiffuseTexture, hitPoint.xz * vec2(0.5,-0.5)).rgb, vec3(2.2));
			}

			bounceIsSpecular = false;
			colorMask *= Tr;
			colorMask *= hitColor;

			// pre-calculate Lambert(cosine-weighted) diffuse lighting amount
			diffuseCosWeight = max(0.0, dot(dirToLight, nl));

			if (hitShininess > 0.0)
			{
				// pre-calculate Blinn-Phong(powered cosine) specular lighting amount
				halfDirection = normalize(-rayDirection + dirToLight);
				specularCosWeight = max(0.0, dot(halfDirection, nl));
				specularCosWeight = pow(specularCosWeight, hitShininess);
			}

			// create shadow ray
			rayDirection = dirToLight;
			rayOrigin = hitPoint + nl * 0.01;

			isShadowRay = true;
			//shadowRayAimedAtPositionalLight = true;
			continue;
			
		} // end if (hitType == CLEARCOAT_DIFFUSE || hitType == CHECKER)

	} // end for (int bounces = 0; bounces < 6; bounces++)

	return max(vec3(0), finalColor);
} // end vec3 CalculateRadiance()


void DefineScene()
{
	directionalLights[0] = DirectionalLight(vec3(-1.0, 1.0,-1.0), vec3(1.0, 1.0, 1.0), 1.0, DIRECTIONAL_LIGHT);

	pointLights[0] = PointLight(vec3(2.0, 5.0, 5.0), vec3(1.0, 1.0, 1.0), 1.0, POINT_LIGHT);

	float angle = mod(uTime, TWO_PI);
	vec3 spotLightPosition = vec3(0.0, 5.0, 0.0);
	vec3 spotLightTarget = vec3(cos(angle) * 5.0, 0.0, sin(angle) * 5.0);
	vec3 spotLightAimDirection = normalize(spotLightTarget - spotLightPosition); 
	spotLights[0] = SpotLight(spotLightPosition, spotLightAimDirection, 0.5, vec3(1.0, 1.0, 1.0), 1.0, POINT_LIGHT);
} // end void DefineScene()



float tentFilter(float x)
{
	return (x < 0.5) ? sqrt(2.0 * x) - 1.0 : 1.0 - sqrt(2.0 - (2.0 * x));
}

void main()
{
	/*
	// simple shader-only camera implementation (commented out)

	vec3 worldUp = vec3(0, 1, 0);
	vec3 camRight, camUp, camForward;
	vec3 lookTarget = vec3(0,0,0);
	vec3 lookDirection = normalize( lookTarget - cameraPosition );
	// camera looks down -Z axis, but the 'camForward' basis vector must point the opposite way, so that it points to +Z axis
	// in the end, we need all camera basis vectors pointing to positive +X, +Y, and +Z axes: 
	// camRight points along the +X axis, camUp points along the +Y axis, 'camForward' points along the +Z axis
	camForward = -lookDirection; // by flipping this around, we get back to the true 'camForward' basis vector, (even though the camera 'looks' in the opposite Z direction)
	camRight = cross(worldUp, camForward);
	camRight = normalize(camRight);
	camUp = cross(camForward, camRight);
	camUp = normalize(camUp); 
	
	*/

	vec3 camRight   = vec3( uCameraMatrix[0][0],  uCameraMatrix[0][1],  uCameraMatrix[0][2]);
	vec3 camUp      = vec3( uCameraMatrix[1][0],  uCameraMatrix[1][1],  uCameraMatrix[1][2]);
	vec3 camForward = vec3(-uCameraMatrix[2][0], -uCameraMatrix[2][1], -uCameraMatrix[2][2]);
	vec3 cameraPos  = vec3( uCameraMatrix[3][0],  uCameraMatrix[3][1],  uCameraMatrix[3][2]);

	// calculate unique seed for rng() function
	seed = uvec2(uFrameCounter, uFrameCounter + 1.0) * uvec2(gl_FragCoord);

	vec2 pixelOffset = vec2(tentFilter(rng()), tentFilter(rng())) * 0.5;
	//pixelOffset = vec2(0);
	vec2 uv = (gl_FragCoord.xy + pixelOffset) / uResolution;
	vec2 pixelPos = uv * 2.0 - 1.0;
	vec3 rayDir = normalize( pixelPos.x * camRight * uULen + pixelPos.y * camUp * uVLen + camForward );

	rayOrigin = cameraPos;
	rayDirection = rayDir;

	DefineScene();

	// calculate pixel color through ray tracing
	vec3 currentPixelColor = CalculateRadiance();

	vec3 previousPixelColor = texelFetch(uPreviousScreenImageTexture, ivec2(gl_FragCoord.xy), 0).rgb;

	// for dynamic scenes
	if (uSceneIsDynamic)
	{
		if (uCameraIsMoving) // camera is currently moving
		{
			previousPixelColor *= 0.5; // motion-blur trail amount (old image)
			currentPixelColor *= 0.5; // brightness of new image (noisy)
		}
		else
		{
			previousPixelColor *= 0.6; // motion-blur trail amount (old image)
			currentPixelColor *= 0.4; // brightness of new image (noisy)
		}
	}

	// for static scenes
	if (!uSceneIsDynamic)
	{
		if (uFrameCounter == 1.0) // camera just moved after being still
		{
			previousPixelColor = vec3(0); // clear rendering accumulation buffer
		}
		else if (uCameraIsMoving) // camera is currently moving
		{
			previousPixelColor *= 0.5; // motion-blur trail amount (old image)
			currentPixelColor *= 0.5; // brightness of new image (noisy)
		}
	}
	
	pixelColor = vec4(previousPixelColor + currentPixelColor, 1);

} // end void main()
`;

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
uniformLocation_uDiffuseTexture = gl.getUniformLocation(rayTracingShaderProgram, 'uDiffuseTexture');
uniformLocation_uMatrices = gl.getUniformLocation(rayTracingShaderProgram, 'uMatrices');
uniformLocation_uCameraMatrix = gl.getUniformLocation(rayTracingShaderProgram, 'uCameraMatrix');
uniformLocation_uSphere0InvMatrix = gl.getUniformLocation(rayTracingShaderProgram, 'uSphere0InvMatrix');
uniformLocation_uBox0InvMatrix = gl.getUniformLocation(rayTracingShaderProgram, 'uBox0InvMatrix');
uniformLocation_uCameraIsMoving = gl.getUniformLocation(rayTracingShaderProgram, 'uCameraIsMoving');
uniformLocation_uSceneIsDynamic = gl.getUniformLocation(rayTracingShaderProgram, 'uSceneIsDynamic');
uniformLocation_uTime = gl.getUniformLocation(rayTracingShaderProgram, 'uTime');
uniformLocation_uFrameCounter = gl.getUniformLocation(rayTracingShaderProgram, 'uFrameCounter');
uniformLocation_uResolution = gl.getUniformLocation(rayTracingShaderProgram, 'uResolution');
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
	gl.uniform1f(uniformLocation_uCameraIsMoving, uCameraIsMoving);
	gl.uniform1f(uniformLocation_uTime, uTime);
	gl.uniform1f(uniformLocation_uFrameCounter, uFrameCounter);

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

	// BOX0
	// the following clears out object's matrix from last frame
	box0.updateMatrixWorld(true);
	// now build up a series of transformations on the object (position, rotation, scale, shear)
	box0.position.set(Math.sin(uTime) * 2, 1, 0);
	//box0.position.set(-3, 1, 0);

	//box0.rotation.set(0, uTime % TWO_PI, 0);

	scale = Math.abs(Math.sin(uTime)) * 2 + 0.1;
	box0.scale.set(1, scale, 1);

	shearMatrix.makeShear(0.0, 0.0, 0.0, 0.0, Math.sin(uTime), 0.0); // (xy, xz,  yx, yz,  zx, zy)
	box0.matrixWorld.multiply(shearMatrix);

	// finally, send box0's Matrix Inverse to the GPU
	uBox0InvMatrix.copy(box0.matrixWorld).invert();
	gl.uniformMatrix4fv(uniformLocation_uBox0InvMatrix, false, uBox0InvMatrix.elements);



	gl.activeTexture(gl.TEXTURE0);
	gl.bindTexture(gl.TEXTURE_2D, uPreviousScreenImageTexture);
	// Tell the shader to use texture unit 0 for uPreviousScreenImageTexture
	gl.uniform1i(uniformLocation_uPreviousScreenImageTexture, 0);
	

	gl.activeTexture(gl.TEXTURE1);
	gl.bindTexture(gl.TEXTURE_2D, uDiffuseTexture);
	// Tell the shader to use texture unit 1 for uDiffuseTexture
	gl.uniform1i(uniformLocation_uDiffuseTexture, 1);

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
