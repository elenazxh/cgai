
/////////////////////////////////////////////////////
//// CS 8803/4803 CGAI: Computer Graphics in AI Era
//// Assignment 2A: Volumetric Ray Tracing
/////////////////////////////////////////////////////

precision highp float;              //// set default precision of float variables to high precision

varying vec2 vUv;                   //// screen uv coordinates (varying, from vertex shader)
uniform vec2 iResolution;           //// screen resolution (uniform, from CPU)
uniform float iTime;                //// time elapsed (uniform, from CPU)
uniform highp sampler3D iVolume;    //// volume texture

/////////////////////////////////////////////////////
//// camera initialization
/////////////////////////////////////////////////////

//// set camera: ro - camera position, ta - camera lookat, cr - camera rotation angle
mat3 setCamera(in vec3 ro, in vec3 ta, float cr)
{
	vec3 cw = normalize(ta-ro);
	vec3 cp = vec3(sin(cr), cos(cr),0.0);
	vec3 cu = normalize(cross(cw,cp));
	vec3 cv = cross(cu,cw);
    return mat3(cu, cv, cw);
}

/////////////////////////////////////////////////////
//// density-to-color conversion
/////////////////////////////////////////////////////

//// Inigo Quilez - https://iquilezles.org/articles/palettes/
//// This function converts a scalar value t to a color
vec3 palette(in float t) 
{
  vec3 a = vec3(0.5, 0.5, 0.5);
  vec3 b = vec3(0.5, 0.5, 0.5);
  vec3 c = vec3(1.0, 1.0, 1.0);
  vec3 d = vec3(0.0, 0.10, 0.20);

  return a + b * cos(6.28318 * (c * t + d));
}

/////////////////////////////////////////////////////
//// sdf definitions
/////////////////////////////////////////////////////

//// sdf sphere
float sdSphere(vec3 p, float s)
{
    return length(p)-s;
}

//// sdf box
float sdBox(vec3 p, vec3 b)
{
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

/////////////////////////////////////////////////////
//// color and density calculation from volume data
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//// Step 1: calculate color and density from sdf
//// You are asked to convert the negative sdf value to a vec4 with the first three components as color rgb values and the last component as density.
//// For color, you may use the provided palette function to convert the sdf value to an rgb color.
//// For density, we assume it is alway 1.0 inside the object (sdf < 0) and 0.0 outside the object (sdf >= 0).
/////////////////////////////////////////////////////

vec4 readSDFVolume(vec3 p)
{
    //// sdf object
    float distance = sdSphere(p, 1.0); 

    //// convert sdf value to a color

    //// your implementation starts

    vec3 color = palette(-distance);
    float density = distance < 0.0 ? 1.0 : 0.0;

    return vec4(color, density);

    //// your implementation ends
}

/////////////////////////////////////////////////////
//// Step 2: calculate color and density from CT data
//// You are asked to convert the CT data to a vec4 with the first three components as color rgb values and the last component as density.
//// For density, you should read the density value from the first component of the volumetric texture iVolume with tex_coord.
//// For color, you should use the provided palette function to convert the density value to an rgb color.
//// You may want to multiple the returned vec4 with a constant to enhance the visualization color.
/////////////////////////////////////////////////////

vec4 readCTVolume(vec3 p)
{
    //// normalize coordinates to [0, 1] range
    vec3 tex_coord = (p + vec3(1.0)) * 0.5;
    //// check if tex_coord is outside the box
    if (tex_coord.x < 0.0 || tex_coord.x > 1.0 || 
        tex_coord.y < 0.0 || tex_coord.y > 1.0 || 
        tex_coord.z < 0.0 || tex_coord.z > 1.0) {
        return vec4(0.0);
    }

    //// your implementation starts

    float density = texture(iVolume, tex_coord).r;
    vec3 color = palette(density);

    return vec4(color, density) * 2.0;

    //// your implementation ends
}

vec4 readMoonVolume(vec3 p)
{
    //// sdf object
    float distance = sdSphere(p, 1.0); 

    //// convert sdf value to a color

    //// your implementation starts

    vec3 color_a = vec3(0.);
    vec3 color_b = vec3(1.);

    vec3 color = mix(vec3(0.984, 0.882, 0.816), vec3(0.), p.y-0.6);

    float density = distance < 0.0 ? 1.0 : 0.0;

    return vec4(color, density);

    //// your implementation ends
}
/*
float sdMoon(vec2 p, float r, float d) {
    return max(length(p) - r, -length(p + vec2(d, 0.0)) + r);
}
vec4 readMoonVolume(vec3 p)
{
    //// sdf object
    float distance = sdMoon(p.xy, 1.5, 0.2);; 

    //// convert sdf value to a color

    //// your implementation starts

    vec3 color = palette(-distance);
    float density = distance < 0.0 ? 1.0 : 0.0;

    return vec4(color, density);

    //// your implementation ends
}
*/
/*
float hash(vec3 p) {
    return fract(sin(dot(p, vec3(127.1, 311.7, 74.7))) * 4375.86713);
}
// Simplex noise
float noise(vec3 p) {
    vec3 i = floor(p);
    vec3 f = fract(p);
    
    vec3 u = f * f * (3.0 - 2.0 * f);
    
    return mix(mix(mix(hash(i + vec3(0, 0, 0)), hash(i + vec3(1, 0, 0)), u.x),
                   mix(hash(i + vec3(0, 1, 0)), hash(i + vec3(1, 1, 0)), u.x), u.y),
               mix(mix(hash(i + vec3(0, 0, 1)), hash(i + vec3(1, 0, 1)), u.x),
                   mix(hash(i + vec3(0, 1, 1)), hash(i + vec3(1, 1, 1)), u.x), u.y), u.z);
}
*/

/* Value noise derivative from Inigo Quilez */
float hash(float p) {
    p = fract(p * .1031);
    p *= p + 33.33;
    p *= p + p;
    return fract(p);
}
float noise( in vec3 x ) {
    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f * f * (3.0 - 2.0 * f);

    float n = p.x + p.y * 57.0 + 113.0 * p.z;

    return mix(mix(mix( hash(n+ 0.0), hash(n+  1.0),f.x),
                        mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y),
                    mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                        mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z);
}

float fbm(vec3 p) {
    float total = 0.0, amplitude = 1.0;
    for (int i = 0; i < 4; i++) {
        total += noise(p) * amplitude;
        p *= 2.0;
        amplitude *= 0.5;
    }
    return total;
}
vec4 readCloudVolume(vec3 p)
{
    float heightFactor = smoothstep(0.0, 1.0, p.y) * smoothstep(1.5, 0.5, p.y);
    float density = fbm(p * 3.0) * heightFactor;
    density = max(density - 0.3, 0.0) * 2.5;

    vec3 color = mix(vec3(1.0), vec3(0.76, 0.86, 1.0), density);

    return vec4(color, density);
}

/////////////////////////////////////////////////////
//// Step 3: ray tracing with volumetric data
//// You are asked to implement the front-to-back volumetric ray tracing algorithm to accummulate colors along each ray. 
//// Your task is to accumulate color and transmittance along the ray based on the absorption-emission volumetric model.
//// You may want to read the course slides, Equation (3) in the original NeRF paper, and the A2a document for the rendering model.
/////////////////////////////////////////////////////

//// ro - ray origin, rd - ray direction, 
//// near - near bound of t, far - far bound of t, 
//// n_samples - number of samples between near and far
vec4 volumeRendering(vec3 ro, vec3 rd, float near, float far, int n_samples) 
{
    float stepSize = (far - near) / float(n_samples);                           //// step size of each sample

    //// color and transmittance
    vec3 color = vec3(0.0);                                                     //// accumulated color
    float transmittance = 1.0;                                                  //// transmittance

    //// ray marching loop
    for (int i = 0; i < n_samples; i++) {
        float t = near + stepSize * float(i);                                   //// t value along the ray
        vec3 p = ro + t * rd;                                                   //// sample position on the ray

        //// your implementation starts

        /* default rendering
        vec4 sdfSample = readSDFVolume(p - vec3(-2.0, 0.0, 0.0));
        vec4 ctSample = readCTVolume(p - vec3(2.0, 0.0, 0.0));

        float sigma = sdfSample.a + ctSample.a;
        vec3 c_i = sdfSample.rgb + ctSample.rgb;
        */ 

        // Cloud Volume Rendering
        vec4 moonSample = readMoonVolume(p - vec3(0., 0.85,0.));
        vec4 cloudSample = readCloudVolume(p);
        float sigma = moonSample.a + cloudSample.a;
        vec3 c_i = (moonSample.rgb * moonSample.a + cloudSample.rgb * cloudSample.a) / max(sigma, 0.0001);

        color += transmittance * (1.0 - exp(-sigma * stepSize)) * c_i;
        transmittance *= exp(-sigma * stepSize);

        //// your implementation ends

        //// early termination if opacity is high
        if (transmittance < 0.01) break;
    }
    
    //// return color and transmittance
    return vec4(color, 1.0 - transmittance);                                   
}


void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    //// normalize fragment coordinates to [-0.5, 0.5] range
    vec2 uv = (fragCoord - 0.5 * iResolution.xy) / iResolution.y;

    //// camera 
    float angle = 0.5 * iTime;           //0.1 for cloud                                   //// camera angle
    vec3 ta = vec3(0.0, 0.0, 0.0);                                              //// object center
    float radius = 5.5;                                                         //// camera rotation
    float height = 2.2;                                                         //// camera height
    vec3 ro = ta + vec3(radius * cos(angle), height, radius * sin(angle));      //// camera position
    mat3 ca = setCamera(ro, ta, 0.0);                                           //// camera matrix
    
    //// ray
    vec3 rd = ca * normalize(vec3(uv, 1.0));                                    //// ray direction
    float near = 2.0;                                                           //// near bound
    float far = 8.5;                                                            //// far bound    
    int n_samples = 256;                                                        //// number of samples along each ray
    
    //// final output color
    fragColor = volumeRendering(ro, rd, near, far, n_samples);                  //// volumetric ray marching
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}