// common Vertex shader source code
const commonVertCode =
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

//screenCopyFragShader source code
const screenCopyFragCode =
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
const screenOutputFragCode =
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

	pixelColor.rgb *= uOneOverSampleCounter;

	// apply tone mapping (brings pixel into 0.0-1.0 rgb color range)
	pixelColor.rgb = ReinhardToneMapping(pixelColor.rgb);

	// lastly, apply gamma correction (gives more intensity/brightness range where it's needed)
	pixelColor = clamp(vec4( pow(pixelColor.rgb, vec3(0.4545)), 1.0 ), 0.0, 1.0);
}
`;


//rayTracingFragShader source code
const rayTracingFragCode =
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
uniform sampler2D uBlueNoiseTexture;
uniform sampler2D uDiffuseTexture;
uniform mat4 uMatrices[2];
uniform mat4 uCameraMatrix;
uniform mat4 uSphere0InvMatrix;
uniform mat4 uSphere1InvMatrix;
uniform mat4 uBox0InvMatrix;
uniform mat4 uBox1InvMatrix;
uniform vec2 uResolution;
uniform vec2 uRandomVec2;
uniform float uTime;
uniform float uFrameCounter;
uniform float uULen;
uniform float uVLen;
uniform float uFocusDistance;
uniform float uApertureSize;
uniform bool uSceneIsDynamic;
uniform bool uCameraIsMoving;

out vec4 pixelColor; // final pixel color output

// global variables
vec3 rayOrigin; 
vec3 rayDirection;
uvec2 seed;

struct DirectionalLight{ vec3 directionTowardsLight; vec3 color; float power; int type; };
struct PointLight{ vec3 position; vec3 color; float power; int type; };
struct SpotLight{ vec3 position; vec3 spotLightAimingDirection; float cutoffAngle; vec3 color; float power; int type; };

DirectionalLight directionalLights[N_DIRECTIONAL_LIGHTS];
PointLight pointLights[N_POINT_LIGHTS];
SpotLight spotLights[N_SPOT_LIGHTS];


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

// globals used in rand() function
vec4 randVec4; // samples and holds the RGBA blueNoise texture value for this pixel
float randNumber; // the final randomly generated number (range: 0.0 to 1.0)
float counter; // will get incremented by 1 on each call to rand()
int channel; // the final selected color channel to use for rand() calc (range: 0 to 3, corresponds to R,G,B, or A)
float rand()
{
	counter++; // increment counter by 1 on every call to rand()
	// cycles through channels, if modulus is 1.0, channel will always be 0 (the R color channel)
	channel = int(mod(counter, 2.0)); 
	// but if modulus was 4.0, channel will cycle through all available channels: 0,1,2,3,0,1,2,3, and so on...
	randNumber = randVec4[channel]; // get value stored in channel 0:R, 1:G, 2:B, or 3:A
	return fract(randNumber); // we're only interested in randNumber's fractional value between 0.0 (inclusive) and 1.0 (non-inclusive)
}

vec3 randomDirectionInPhongSpecular(vec3 reflectionDir, float shininess)
{
	shininess *= 50.0;
	float cosTheta = pow(rng(), 1.0 / (2.0 + shininess));
	float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
	float phi = rng() * TWO_PI;
	float x = sinTheta * cos(phi);
	float y = sinTheta * sin(phi);

	vec3 U = normalize( cross( abs(reflectionDir.y) < 0.9 ? vec3(0, 1, 0) : vec3(0, 0, 1), reflectionDir ) );
	vec3 V = cross(reflectionDir, U);
	return normalize(x * U + y * V + cosTheta * reflectionDir);
}


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
		hitColor = vec3(1.0, 1.0, 1.0);
		hitShininess = 1000.0;
		hitType = CHECKER;
	}

	// transform ray into Ellipsoid Param's object space
	rObjOrigin = vec3( uSphere0InvMatrix * vec4(rayOrigin, 1.0) );
	rObjDirection = vec3( uSphere0InvMatrix * vec4(rayDirection, 0.0) );
	
	float angleAmount = (sin(uTime) * 0.5 + 0.5);
	d = EllipsoidParamIntersect(-1.0, 1.0, TWO_PI, rObjOrigin, rObjDirection, n);
	//d = EllipsoidParamIntersect(-1.0, angleAmount, PI, rObjOrigin, rObjDirection, n);
	if (d < t)
	{
		t = d;
		hitNormal = (transpose(mat3(uSphere0InvMatrix)) * n);
		hitColor = vec3(1);//vec3(0.0, 0.3, 1.0);
		hitShininess = 1000.0;
		hitType = TRANSPARENT;
	}

	// transform ray into Ellipsoid Param's object space
	rObjOrigin = vec3( uSphere1InvMatrix * vec4(rayOrigin, 1.0) );
	rObjDirection = vec3( uSphere1InvMatrix * vec4(rayDirection, 0.0) );
	
	d = EllipsoidParamIntersect(-1.0, 1.0, TWO_PI, rObjOrigin, rObjDirection, n);
	if (d < t)
	{
		t = d;
		hitNormal = (transpose(mat3(uSphere1InvMatrix)) * n);
		hitColor = vec3(1,0.6,0.1);
		hitShininess = 100.0;
		hitType = METAL;
	}

	// transform ray into Unit Box0's object space
	rObjOrigin = vec3( uBox0InvMatrix * vec4(rayOrigin, 1.0) );
	rObjDirection = vec3( uBox0InvMatrix * vec4(rayDirection, 0.0) );
	
	d = UnitBoxIntersect(rObjOrigin, rObjDirection, n);
	if (d < t)
	{
		t = d;
		hitNormal = (transpose(mat3(uBox0InvMatrix)) * n);
		hitColor = vec3(0.0, 1.0, 0.0);
		hitShininess = 1000.0;
		hitType = DIFFUSE;
	}

	// transform ray into Unit Box1's object space
	rObjOrigin = vec3( uBox1InvMatrix * vec4(rayOrigin, 1.0) );
	rObjDirection = vec3( uBox1InvMatrix * vec4(rayDirection, 0.0) );
	
	d = UnitBoxIntersect(rObjOrigin, rObjDirection, n);
	if (d < t)
	{
		t = d;
		hitNormal = (transpose(mat3(uBox1InvMatrix)) * n);
		hitColor = vec3(1.0, 1.0, 0.0);
		hitShininess = 1000.0;
		hitType = DIFFUSE;
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

	/*
	float angleBetweenRayAndLight = acos(dot(rayDirection, -spotLights[0].spotLightAimingDirection));
	if (angleBetweenRayAndLight > spotLights[0].cutoffAngle)
		return t;

	d = BoundingBoxIntersect(spotLights[0].position - vec3(0.1, 0.1, 0.25), spotLights[0].position + vec3(0.1, 0.1, 0.25), rayOrigin, 1.0/rayDirection);
	if (d < t)
	{
		t = d;
		hitColor = spotLights[0].color * spotLights[0].power;
		hitType = SPOT_LIGHT;
	}
 	*/
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
	int hitType = -100;

	vec3 checkColor0 = vec3(0.1);
	vec3 checkColor1 = vec3(1.0);
	vec3 colorMask = vec3(1);
	vec3 finalColor = vec3(0);
	vec3 ambientColor = vec3(0);
	vec3 diffuseColor = vec3(0);
	vec3 specularColor = vec3(0);
	vec3 n, nl, hitPoint;
	vec3 dirToLight;
	vec3 halfDirection;
	vec3 tdir;

	float diffuseCosWeight = 0.0;
	float specularCosWeight = 0.0;
	float checkScale = 1.0;
	float t = INFINITY;
	float nc, nt, ratioIoR, Re, Tr;
	float P, RP, TP;
	float distanceSquaredToLight;

	bool bounceIsSpecular = true;
	bool isShadowRay = false;
	bool isCausticRay = false;
	bool shadowRayAimedAtPositionalLight = false;


	for (int bounces = 0; bounces < 5; bounces++)
	{
		t = SceneIntersect(rayOrigin, rayDirection, shadowRayAimedAtPositionalLight, 
					hitNormal, hitColor, hitShininess, hitType);
		
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

				finalColor = diffuseColor + specularColor; // diffuse color already includes ambient color in this method

				if (isCausticRay)
					finalColor = diffuseColor * pow(max(dot(-normalize(hitNormal), normalize(dirToLight)), 0.0), 1.5) * 1.5;
			}

			break;
		} // end if (t == INFINITY && !shadowRayAimedAtPositionalLight)


		if (hitType == POINT_LIGHT || hitType == SPOT_LIGHT)
		{
			// set object's ambient lighting (this is the darkest it can be, even in shadow)
			ambientColor = colorMask * 0.2;

			diffuseColor = colorMask * hitColor;
			diffuseColor = mix(ambientColor, diffuseColor, diffuseCosWeight);

			specularColor = specularCosWeight * hitColor;

			finalColor = diffuseColor + specularColor; // diffuse color already includes ambient color in this method

			break;
		} // end if (hitType == POINT_LIGHT || hitType == SPOT_LIGHT)

		
		// if still is shadow ray, the ray aimed towards the light is blocked by a scene object (it hit something), 
		// and therefore the surface remains in shadow 
		if (hitType != TRANSPARENT && isShadowRay)
		{
			// ambientColor is for darker parts of the scene that are unlit by direct light sources
			// set object's ambient lighting (this is the darkest it can be, even in shadow)
			ambientColor = colorMask * 0.2;

			finalColor = ambientColor;

			break;
		}
		if (hitType == TRANSPARENT && isShadowRay)
			isCausticRay = true;
			

		// useful data
		n = normalize(hitNormal);
		nl = dot(n, rayDirection) < 0.0 ? n : -n;
		hitPoint = rayOrigin + rayDirection * t;
		
		//dirToLight = directionalLights[0].directionTowardsLight;
		dirToLight = normalize(pointLights[0].position - hitPoint);
		/* 
		dirToLight = spotLights[0].position - hitPoint;
		distanceSquaredToLight = max(1.0, dot(dirToLight, dirToLight) * 0.2);
		dirToLight = normalize(dirToLight);
 		*/

		if (hitType == DIFFUSE)
		{
			bounceIsSpecular = false;
			colorMask *= hitColor;

			// pre-calculate Lambert(cosine-weighted) diffuse lighting amount
			diffuseCosWeight = max(0.0, dot(dirToLight, nl));
			///diffuseCosWeight /= distanceSquaredToLight;
			
			if (hitShininess > 0.0)
			{
				// pre-calculate Blinn-Phong(powered cosine) specular lighting amount
				halfDirection = normalize(-rayDirection + dirToLight);
				specularCosWeight = max(0.0, dot(halfDirection, nl));
				specularCosWeight = pow(specularCosWeight, hitShininess);
				///specularCosWeight /= distanceSquaredToLight;
			}
			
			// create shadow ray
			rayDirection = dirToLight;
			rayOrigin = hitPoint + nl * 0.01;

			isShadowRay = true;
			shadowRayAimedAtPositionalLight = true;
			continue;
		} // end if (hitType == DIFFUSE)


		if (hitType == METAL)
		{
			colorMask *= hitColor;

			if (hitShininess > 0.0)
			{
				// pre-calculate Blinn-Phong(powered cosine) specular lighting amount
				halfDirection = normalize(-rayDirection + dirToLight);
				specularCosWeight = max(0.0, dot(halfDirection, nl));
				specularCosWeight = pow(specularCosWeight, hitShininess);
			}
			
			// create reflection ray
			vec3 reflectedRay = reflect(rayDirection, nl);
			//vec3 randomizedRay = normalize( reflectedRay + vec3(rand() * 2.0 - 1.0, rand() * 2.0 - 1.0, rand() * 2.0 - 1.0) );
			vec3 randomizedRay = randomDirectionInPhongSpecular(reflectedRay, hitShininess);
			rayDirection = randomizedRay;
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

			if (bounces == 0 && rand() < P)
			{
				if (hitShininess > 0.0)
				{
					// pre-calculate Blinn-Phong(powered cosine) specular lighting amount
					halfDirection = normalize(-rayDirection + dirToLight);
					specularCosWeight = max(0.0, dot(halfDirection, nl));
					specularCosWeight = pow(specularCosWeight, hitShininess);
				}

				colorMask *= RP;
				rayDirection = reflect(rayDirection, nl); // reflect ray from surface
				rayOrigin = hitPoint + nl * 0.01;
				continue;
			}

			// transmit ray through surface
			colorMask *= hitColor;
			colorMask *= TP;
			
			tdir = refract(rayDirection, nl, ratioIoR);
			if (isCausticRay)
				tdir = rayDirection;
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

			if (rand() < P)
			{
				if (hitShininess > 0.0)
				{
					// pre-calculate Blinn-Phong(powered cosine) specular lighting amount
					halfDirection = normalize(-rayDirection + dirToLight);
					specularCosWeight = max(0.0, dot(halfDirection, nl));
					specularCosWeight = pow(specularCosWeight, hitShininess);
					///specularCosWeight /= distanceSquaredToLight;
				}

				colorMask *= RP;
				rayDirection = reflect(rayDirection, nl); // reflect ray from surface
				rayOrigin = hitPoint + nl * 0.01;
				continue;
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
			colorMask *= TP;
			colorMask *= hitColor;

			// pre-calculate Lambert(cosine-weighted) diffuse lighting amount
			diffuseCosWeight = max(0.0, dot(dirToLight, nl));
			///diffuseCosWeight /= distanceSquaredToLight;

			if (hitShininess > 0.0)
			{
				// pre-calculate Blinn-Phong(powered cosine) specular lighting amount
				halfDirection = normalize(-rayDirection + dirToLight);
				specularCosWeight = max(0.0, dot(halfDirection, nl));
				specularCosWeight = pow(specularCosWeight, hitShininess);
				///specularCosWeight /= distanceSquaredToLight;
			}

			// create shadow ray
			rayDirection = dirToLight;
			rayOrigin = hitPoint + nl * 0.01;

			isShadowRay = true;
			shadowRayAimedAtPositionalLight = true;
			continue;
			
		} // end if (hitType == CLEARCOAT_DIFFUSE || hitType == CHECKER)

	} // end for (int bounces = 0; bounces < 5; bounces++)

	return max(vec3(0), finalColor);
} // end vec3 CalculateRadiance()


void DefineScene()
{
	float angle = mod(uTime, TWO_PI);

	directionalLights[0] = DirectionalLight(normalize(vec3(cos(angle), 1.0, sin(angle))), vec3(1.0, 1.0, 1.0), 2.0, DIRECTIONAL_LIGHT);

	pointLights[0] = PointLight(vec3(cos(angle) * 5.0, 5.0, sin(angle) * 5.0), vec3(1.0, 1.0, 1.0), 1.0, POINT_LIGHT);

	/*
	vec3 spotLightPosition = vec3(0, 10, 0);
	vec3 spotLightTarget = vec3(cos(angle) * 5.0, 0.0, sin(angle) * 5.0);
	vec3 spotLightAimDirection = normalize(spotLightTarget - spotLightPosition); 
	spotLights[0] = SpotLight(spotLightPosition, spotLightAimDirection, 0.3, vec3(1.0, 1.0, 1.0), 50.0, SPOT_LIGHT); 
	*/
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

	// initialize rand() variables
	counter = -1.0; // will get incremented by 1 on each call to rand()
	channel = 0; // the final selected color channel to use for rand() calc (range: 0 to 3, corresponds to R,G,B, or A)
	randNumber = 0.0; // the final randomly-generated number (range: 0.0 to 1.0)
	randVec4 = vec4(0); // samples and holds the RGBA blueNoise texture value for this pixel
	randVec4 = texelFetch(uBlueNoiseTexture, ivec2(mod(gl_FragCoord.xy + floor(uRandomVec2 * 256.0), 256.0)), 0);

	// calculate unique seed for rng() function
	seed = uvec2(uFrameCounter, uFrameCounter + 1.0) * uvec2(gl_FragCoord);

	vec2 pixelOffset = vec2(tentFilter(rand()), tentFilter(rand())) * 0.5;
	vec2 uv = (gl_FragCoord.xy + pixelOffset) / uResolution;
	vec2 pixelPos = uv * 2.0 - 1.0;
	vec3 rayDir = normalize( pixelPos.x * camRight * uULen + pixelPos.y * camUp * uVLen + camForward );

	// depth of field
	vec3 focalPoint = uFocusDistance * rayDir;
	float randomAngle = rand() * TWO_PI; // pick random point on aperture
	float randomRadius = rand() * uApertureSize;
	vec3  randomAperturePos = ( cos(randomAngle) * camRight + sin(randomAngle) * camUp ) * sqrt(randomRadius);
	// point on aperture to focal point
	vec3 finalRayDir = normalize(focalPoint - randomAperturePos);

	rayOrigin = cameraPos + randomAperturePos;
	rayDirection = finalRayDir;


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
			previousPixelColor *= 0.7; // motion-blur trail amount (old image)
			currentPixelColor *= 0.3; // brightness of new image (noisy)
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