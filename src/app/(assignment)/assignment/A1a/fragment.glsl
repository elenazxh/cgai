/////////////////////////////////////////////////////
//// CS 8803/4803 CGAI: Computer Graphics in AI Era
//// Assignment 1A: SDF and Ray Marching
/////////////////////////////////////////////////////

precision highp float;              //// set default precision of float variables to high precision

varying vec2 vUv;                   //// screen uv coordinates (varying, from vertex shader)
uniform vec2 iResolution;           //// screen resolution (uniform, from CPU)
uniform float iTime;                //// time elapsed (uniform, from CPU)

const vec3 CAM_POS = vec3(-0.35, 1.2, -3.0);

/////////////////////////////////////////////////////
//// sdf functions
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//// Step 1: sdf primitives
//// You are asked to implement sdf primitive functions for sphere, plane, and box.
//// In each function, you will calculate the sdf value based on the function arguments.
/////////////////////////////////////////////////////

//// sphere: p - query point; c - sphere center; r - sphere radius
float sdfSphere(vec3 p, vec3 c, float r)
{
    //// your implementation starts
    
    return length(p-c) - r;
    
    //// your implementation ends
}

//// plane: p - query point; h - height
float sdfPlane(vec3 p, float h)
{
    //// your implementation starts
    
    return p.y - h;
    
    //// your implementation ends
}

//// box: p - query point; c - box center; b - box half size (i.e., the box size is (2*b.x, 2*b.y, 2*b.z))
float sdfBox(vec3 p, vec3 c, vec3 b)
{
    //// your implementation starts
    
    vec3 d = abs(p-c) - b;
    return length(max(d, 0.0)) + min(max(d.x, max(d.y, d.z)), 0.0);
    
    //// your implementation ends
}

/////////////////////////////////////////////////////
//// boolean operations
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//// Step 2: sdf boolean operations
//// You are asked to implement sdf boolean operations for intersection, union, and subtraction.
/////////////////////////////////////////////////////

float sdfIntersection(float s1, float s2)
{
    //// your implementation starts
    
    return max(s1, s2);

    //// your implementation ends
}

float sdfUnion(float s1, float s2)
{
    //// your implementation starts
    
    return min(s1, s2);

    //// your implementation ends
}

float sdfSubtraction(float s1, float s2)
{
    //// your implementation starts
    
    return max(s1, -s2);

    //// your implementation ends
}

/////////////////////////////////////////////////////
//// sdf calculation
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//// Step 3: scene sdf
//// You are asked to use the implemented sdf boolean operations to draw the following objects in the scene by calculating their CSG operations.
/////////////////////////////////////////////////////

//// sdf: p - query point
float sdf(vec3 p)
{
    float s = 0.;

    //// 1st object: plane
    float plane1_h = -0.1;
    
    //// 2nd object: sphere
    vec3 sphere1_c = vec3(-2.0, 1.0, 0.0);
    float sphere1_r = 0.25;

    //// 3rd object: box
    vec3 box1_c = vec3(-1.0, 1.0, 0.0);
    vec3 box1_b = vec3(0.2, 0.2, 0.2);

    //// 4th object: box-sphere subtraction
    vec3 box2_c = vec3(0.0, 1.0, 0.0);
    vec3 box2_b = vec3(0.3, 0.3, 0.3);

    vec3 sphere2_c = vec3(0.0, 1.0, 0.0);
    float sphere2_r = 0.4;

    //// 5th object: sphere-sphere intersection
    vec3 sphere3_c = vec3(1.0, 1.0, 0.0);
    float sphere3_r = 0.4;

    vec3 sphere4_c = vec3(1.3, 1.0, 0.0);
    float sphere4_r = 0.3;

    //// calculate the sdf based on all objects in the scene
    
    //// your implementation starts
    
    float s1 = sdfPlane(p, plane1_h);
    float s2 = sdfSphere(p, sphere1_c, sphere1_r);
    float s3 = sdfBox(p, box1_c, box1_b);
    float s4 = sdfSubtraction(sdfBox(p, box2_c, box2_b), sdfSphere(p, sphere2_c, sphere2_r));
    float s5 = sdfIntersection(sdfSphere(p, sphere3_c, sphere3_r), sdfSphere(p, sphere4_c, sphere4_r));

    s = min(s1, min(s2, min(s3, min(s4, s5))));
    
    //// your implementation ends

    return s;
}

/////////////////////////////////////////////////////
//// ray marching
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//// Step 4: ray marching
//// You are asked to implement the ray marching algorithm within the following for-loop.
/////////////////////////////////////////////////////

//// ray marching: origin - ray origin; dir - ray direction 
float rayMarching(vec3 origin, vec3 dir)
{
    float s = 0.0;
    for(int i = 0; i < 100; i++)
    {
        //// your implementation starts

        vec3 p = origin + s * dir;
        s += sdf(p);
        if (sdf(p) < 0.001 || s > 100.0 ) break;

        //// your implementation ends
    }
    
    return s;
}

/////////////////////////////////////////////////////
//// normal calculation
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//// Step 5: normal calculation
//// You are asked to calculate the sdf normal based on finite difference.
/////////////////////////////////////////////////////

//// normal: p - query point
vec3 normal(vec3 p)
{
    float s = sdf(p);          //// sdf value in p
    float dx = 0.01;           //// step size for finite difference

    //// your implementation starts
    
    float delta_x = sdf(p + vec3(dx, 0.0, 0.0)) - sdf(p - vec3(dx, 0.0, 0.0));
    float delta_y = sdf(p + vec3(0.0, dx, 0.0)) - sdf(p - vec3(0.0, dx, 0.0));
    float delta_z = sdf(p + vec3(0.0, 0.0, dx)) - sdf(p - vec3(0.0, 0.0, dx));
    return normalize(vec3(delta_x, delta_y, delta_z));

    //// your implementation ends
}

/////////////////////////////////////////////////////
//// Phong shading
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//// Step 6: lighting and coloring
//// You are asked to specify the color for each object in the scene.
//// Each object must have a separate color without mixing.
//// Notice that we have implemented the default Phong shading model for you.
/////////////////////////////////////////////////////

vec3 phong_shading(vec3 p, vec3 n)
{
    //// background
    if(p.z > 10.0){
        return vec3(0.9, 0.6, 0.2);
    }

    //// phong shading
    vec3 lightPos = vec3(4.*sin(iTime), 4., 4.*cos(iTime));  
    vec3 l = normalize(lightPos - p);               
    float amb = 0.1;
    float dif = max(dot(n, l), 0.) * 0.7;
    vec3 eye = CAM_POS;
    float spec = pow(max(dot(reflect(-l, n), normalize(eye - p)), 0.0), 128.0) * 0.9;

    vec3 sunDir = vec3(0, 1, -1);
    float sunDif = max(dot(n, sunDir), 0.) * 0.2;

    //// shadow
    float s = rayMarching(p + n * 0.02, l);
    if(s < length(lightPos - p)) dif *= .2;

    vec3 color = vec3(1.0, 1.0, 1.0);

    //// your implementation for coloring starts

    if (p.y < 0.1) {
        color = vec3(1.0, 1.0, 0.0);
    } else {
        if (p.x < -1.5) {
            color = vec3(1.0, 0.0, 0.0);
        } 
        else if (p.x < -0.5) {
            color = vec3(0.0, 1.0, 0.0);
        } 
        else if (p.x < 0.5) {
            color = vec3(0.0, 0.0, 1.0);
        } 
        else {
            color = vec3(0.0, 1.0, 1.0);
        }
    }

    //// your implementation for coloring ends

    return (amb + dif + spec + sunDif) * color;
}

/////////////////////////////////////////////////////
//// Step 7: creative expression
//// You will create your customized sdf scene with new primitives and CSG operations in the sdf2 function.
//// Call sdf2 in your ray marching function to render your customized scene.
/////////////////////////////////////////////////////

float smoothUnion(float d1, float d2, float k) {
    float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return mix(d2, d1, h) - k * h * (1.0 - h);
}

//// sdf2: p - query point
// **credit to Inigo Quilez**
float sdRoundedCylinder( vec3 p, float ra, float rb, float h )
{
  vec2 d = vec2( length(p.xz)-2.0*ra+rb, abs(p.y) - h );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0)) - rb;
}
float sdVerticalCapsule( vec3 p, float h, float r )
{
  p.y -= clamp( p.y, 0.0, h );
  return length( p ) - r;
}
float sdCappedCone( vec3 p, float h, float r1, float r2 )
{
  vec2 q = vec2( length(p.xz), p.y );
  vec2 k1 = vec2(r2,h);
  vec2 k2 = vec2(r2-r1,2.0*h);
  vec2 ca = vec2(q.x-min(q.x,(q.y<0.0)?r1:r2), abs(q.y)-h);
  vec2 cb = q - k1 + k2*clamp( dot(k1-q,k2)/dot(k2, k2), 0.0, 1.0 );
  float s = (cb.x<0.0 && ca.y<0.0) ? -1.0 : 1.0;
  return s*sqrt( min(dot(ca, ca),dot(cb, cb)) );
}

float bubbleMotion(vec3 p, float seed)
{
    return sin(iTime * 0.5 * seed + 5.*seed) * 0.5;
}

float sdfBubbles(vec3 p)
{
    float bubble1 = sdfSphere(p, vec3( 0.1, bubbleMotion(p, 2.5) + 1.4, 0), 0.03);
    float bubble2 = sdfSphere(p, vec3(-0.05, bubbleMotion(p, 2.0) + 2.0, 0), 0.035);
    float bubble3 = sdfSphere(p, vec3( 0.15, bubbleMotion(p, 3.0) + 1.8, 0), 0.035);
    float bubble4 = sdfSphere(p, vec3(-0.15, bubbleMotion(p, 3.5) + 1.8, 0), 0.04);

    // Combine all bubbles smoothly
    float bubbles = smoothUnion(bubble1, bubble2, 0.02);
    bubbles = smoothUnion(bubbles, bubble3, 0.02);
    bubbles = smoothUnion(bubbles, bubble4, 0.02);
    
    return bubbles;
}
float sdf_champagne(vec3 p)
{
    float s = 0.;

    float glassBase = sdRoundedCylinder(p - vec3(0.0, -0.1, 0.0), 0.15, 0.02, 0.01);
    float glassStem = sdVerticalCapsule(p - vec3(0.0, -0.1, 0.0), 0.8, 0.04);
    float glassBody1 = sdCappedCone(p - vec3(0.0, 1.555, 0.0), 0.23, 0.3, 0.2);
    float glassBody2 = sdCappedCone(p - vec3(0.0, 1.0, 0.0), 0.32, 0.08, 0.3);
    //float cork = sdTorus(p - vec3(0.0, 1.38, 0.0), vec2(0.2, 0.04));
    

    s = smoothUnion(glassStem, glassBase, 0.1);
    s = smoothUnion(s, glassBody2, 0.1);
    s = smoothUnion(glassBody1, s, 0.01);
    // s = smoothUnion(s, cork, 0.2);

    float bubbles = sdfBubbles(p);

    return smoothUnion(s, bubbles, 0.02);
}

float bubbleMarching(vec3 origin, vec3 dir)
{
    float s = 0.0;
    for(int i = 0; i < 100; i++)
    {
        //// your implementation starts

        vec3 p = origin + s * dir;
        s += sdf_champagne(p);
        if (sdf_champagne(p) < 0.001 || s > 100.0 ) break;

        //// your implementation ends
    }
    
    return s;
}
vec3 bubble_normal(vec3 p)
{
    float s = 0.;          //// sdf value in p
    float dx = 0.01;           //// step size for finite difference

    //// your implementation starts
    
    float delta_x = sdf_champagne(p + vec3(dx, 0.0, 0.0)) - sdf_champagne(p - vec3(dx, 0.0, 0.0));
    float delta_y = sdf_champagne(p + vec3(0.0, dx, 0.0)) - sdf_champagne(p - vec3(0.0, dx, 0.0));
    float delta_z = sdf_champagne(p + vec3(0.0, 0.0, dx)) - sdf_champagne(p - vec3(0.0, 0.0, dx));
    return normalize(vec3(delta_x, delta_y, delta_z));

    //// your implementation ends
}
vec3 bubble_shading(vec3 p, vec3 n)
{
    //// background
    if(p.z > 10.0){
        vec3 color_A = vec3(0.319,0.167,0.365);
        vec3 color_B = vec3(0.628,0.526,0.775);
	    float pct = abs(sin(iTime * 0.7));
    
        vec3 color = mix(color_A, color_B, pct);
        return color;
        // return vec3(.98, .92, .84);
    }

    //// phong shading
    vec3 lightPos = vec3(4.*sin(iTime), 4., 4.*cos(iTime));  
    vec3 l = normalize(lightPos - p);               
    float amb = 0.4;
    float dif = max(dot(n, l), 0.) * 0.7;
    vec3 eye = CAM_POS;

    float shininess = mix(16.0, 128.0, smoothstep(-0.1, 1.7, p.y)); 
    float spec = pow(max(dot(reflect(-l, n), normalize(eye - p)), 0.0), shininess) * 0.9;

    vec3 sunDir = vec3(0, 1, -1);
    float sunDif = max(dot(n, sunDir), 0.) * 0.2;

    //// shadow
    float s = rayMarching(p + n * 0.02, l);
    if(s < length(lightPos - p)) dif *= .2;

    vec3 color = vec3(1.0, 1.0, 1.0);

    //// your implementation for coloring starts

    float bubbleDist = sdfBubbles(p);
    bool isBubble = (bubbleDist < 0.05);

    // Color selection
    if (isBubble) {
        color = vec3(1.0, 0.6, 0.8); // Pink bubbles
    } else {
        float gradientFactor = smoothstep(-0.1, 1.8, p.y);
        color = mix(vec3(0.1, 0.1, 0.0), vec3(0.898, 0.988, 0.988), gradientFactor);
        if (p.y >= 1.78) color = vec3(.98, .92, .84);
    }
    

    //// your implementation for coloring ends

    return (amb + dif + spec + sunDif) * color;
}

/////////////////////////////////////////////////////
//// main function
/////////////////////////////////////////////////////

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - .5 * iResolution.xy) / iResolution.y;         //// screen uv
    vec3 origin = CAM_POS;                                                  //// camera position 
    vec3 dir = normalize(vec3(uv.x, uv.y, 1));                              //// camera direction
    float s = rayMarching(origin, dir);                                     //// ray marching
    vec3 p = origin + dir * s;                                              //// ray-sdf intersection
    vec3 n = normal(p);                                                     //// sdf normal
    vec3 color = phong_shading(p, n);                                       //// phong shading

    fragColor = vec4(color, 1.);                                            //// fragment color
}

void main() 
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}