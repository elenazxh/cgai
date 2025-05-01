// ----------------------------------------------------------------------------------
// Constants and Uniforms
// ----------------------------------------------------------------------------------
precision highp float;              //// set default precision of float variables to high precision

varying vec2 vUv;                   //// screen uv coordinates (varying, from vertex shader)
uniform vec2 iResolution;           //// screen resolution (uniform, from CPU)
uniform float iTime;                //// time elapsed (uniform, from CPU)
// uniform vec4 iMouse;

const vec3 CAM_POS = vec3(-0.35, 1.0, -3.0); //// camera position

struct HitID {
    float dist;
    int id;
};
HitID hit_id = HitID(2000.0, -1);

struct Light {
    vec3 pos;
    vec3 color;
    float ambientStrength;
    float diffuseStrength;
    float specularStrength;
    float shininess;
};

float jitter; // for cloud rendering

float time_cycle; // for time-based animation

// ----------------------------------------------------------------------------------
// Utility Functions
// ----------------------------------------------------------------------------------
vec3 rotateXYZ(vec3 p, vec3 c, vec3 angles)
{   
    p -= c; // Translate to origin
    angles = radians(angles);
    float c1 = cos(angles.x), s1 = sin(angles.x);
    float c2 = cos(angles.y), s2 = sin(angles.y);
    float c3 = cos(angles.z), s3 = sin(angles.z);

    mat3 m = mat3(
        c1 * c3 + s1 * s2 * s3, c2 * s3, -s1 * c3 + c1 * s2 * s3,
        -c1 * s3 + s1 * s2 * c3, c2 * c3, s1 * s3 + c1 * s2 * c3,
        s1 * c2, -s2, c1 * c2
    );

    return m * p + c; // Rotate and translate back
}

float hash(float n)
{
    return fract(sin(n) * 43758.5453);
}

float noise(vec3 x)
{
    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f * f * (3.0 - 2.0 * f);

    float n = p.x + p.y * 57.0 + 113.0 * p.z;

    float res = mix(mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
                        mix(hash(n + 57.0), hash(n + 58.0), f.x), f.y),
                    mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
                        mix(hash(n + 170.0), hash(n + 171.0), f.x), f.y), f.z);
    return res;
}

float fbm(in vec3 p)
{
    mat3 m = mat3(
        0.00,  0.80,  0.60,
    -0.80,  0.36, -0.48,
    -0.60, -0.48,  0.64
    );
    float f;
    f  = 0.5000 * noise(p);  p = m * p * 2.02;
    f += 0.2500 * noise(p);  p = m * p * 2.03;
    f += 0.12500 * noise(p); p = m * p * 2.01;
    f += 0.06250 * noise(p);
    return f;
}

vec3 hash33(vec3 p)
{
    p = fract(p * vec3(443.8975,397.2973, 491.1871));
    p += dot(p.zxy, p.yxz+19.27);
    return fract(vec3(p.x * p.y, p.z*p.x, p.y*p.z));
}

/**
 * Return a float value between 0.0(0s) and 1.0(10s) that represents a 10s cycle
 * Used for time-based animations
 */
float getSyncedTimeCycle() 
{
    float t = float(int(iTime * 120.0) % 1200) / 1200.0; // 10s cycle

    // apply a curve: slow -> fast -> slow
    if (t < 0.5) {
        t = pow(t, 2.0) * 2.0; // Slow in the beginning
    } else {
        t = 1.0 - pow(1.0 - t, 2.0) * 2.0; // Slow in the end
    }
    if (int(iTime) % 20 >= 10) {
        t = 1.0 - t; // Reverse the curve
    }
    return t;
}


// float getSyncedTimeCycle() 
// {
//     // return 0.0;
//     // get mouse position
//     float mouse_x = iMouse.x / iResolution.x;



//     return mouse_x;
// }

// float getSyncedTimeCycle() 
// {
//     // return 0.7;
//     return sin(iTime * 0.2) * 0.5 + 0.5; // 0.0 to 1.0
// }

// ----------------------------------------------------------------------------------
// SDF Shape Functions
// ----------------------------------------------------------------------------------
float sdfSphere(vec3 p, vec3 c, float r)
{
    return length(p - c) - r;
}

float sdfPlane(vec3 p, float h)
{
    return p.y - h;
}

float sdfBox(vec3 p, vec3 c, vec3 b)
{
    vec3 d = abs(p - c) - b;
    return length(max(d, 0.0)) + min(max(d.x, max(d.y, d.z)), 0.0);
}

float sdfEllipsoid(vec3 p, vec3 c, vec3 r)
{
    p = p - c;
    float k0 = length(p / r);
    float k1 = length(p / (r * r));
    return k0 * (k0 - 1.0) / k1;
}

float sdTorus(vec3 p, vec2 t)
{
    vec2 q = vec2(length(p.xz) - t.x, p.y);
    return length(q) - t.y;
}

float sdfOctahedron(vec3 p, vec3 c, float r)
{
    vec3 d = abs(p - c);
    return (d.x + d.y + d.z - r) * 0.57735; // 1/sqrt(3)
}

float sdfOctahedron(in vec3 p, vec3 c, float r, vec3 rot)
{
    p = rotateXYZ(p, c, rot);
    return sdfOctahedron(p, c, r);
}

float sdfVerticleCapsule(in vec3 p, vec3 c, float r, float h)
{
    p -= c;
    p.y -= clamp(p.y, 0.0, h);
    return length(p) - r;
}

float sdfCone(in vec3 p, vec3 c, float angle, float h)
{
    p -= c;
    float q = length(p.xz);
    vec2 ang = vec2(cos(angle), sin(angle));   
    // vec2 ang = vec2(0.980, 0.199);
    return max(dot(ang,vec2(q,p.y)),-h-p.y);
}

// ------------------------------------
// SDF Operations
// ------------------------------------
float sdfIntersection(float s1, float s2)
{
    return max(s1, s2);
}

float sdfUnion(float s1, float s2)
{
    return min(s1, s2);
}

float sdfSubtraction(float s1, float s2)
{
    return max(s1, -s2);
}

float sdfUnionSmooth(float s1, float s2, float k)
{
    return -k * log(exp(-s1 / k) + exp(-s2 / k));
}

float sdfSubtractionSmooth(float s1, float s2, float k)
{
    return -sdfUnionSmooth(-s1, s2, k);
}

// ----------------------------------------------------------------------------------
// Scene-Specific SDFs
// ----------------------------------------------------------------------------------

float sdfCurvyGround(vec3 p, float h) 
{
    p -= vec3(0.0, 0.0, 0.0);
    // float wave = 0.3 * sin(0.5 * p.x) * cos(1.0 * p.z); // Hills and Valleys

    // float texture = 0.02 * sin(40.0 * p.x) * sin(80.0 * p.z); // Texture
    float noise = 0.1 * fbm(p * 20.0) - 0.0;
    // float noise2 = 5.0 * fbm(p * 0.2 - 0.5) - 3.0;
    // Create a height map based on p.y coordinate
    float height = h + abs(p.x) * 0.02;
    return p.y - clamp(height, -0.5, 15.0) - noise;
}

float sdfBorb(vec3 p, vec3 c, float angle, bool birb_hair)
{
    // Rotate around the Y-axis
    float birb_move1 = 3.0 * sin(iTime * 3.0);
    float birb_move2 = 3.0 * cos(iTime * 1.0);
    vec3 local_p = rotateXYZ(p, c, vec3(angle, birb_move1, birb_move2));

    // Define the borb components
    float r = 0.2;
    vec3 head_c  = c + r * vec3(-0.1, 1.5, 0.0);
    vec3 body_c  = c + r * vec3(0.0, 0.0, 0.0);
    vec3 tail_c  = c + r * vec3(1.5, -0.3, 0.0);
    vec3 wingL_c = c + r * vec3(0.0, 0.1, -1.0);
    vec3 wingR_c = c + r * vec3(0.0, 0.1, 1.0);
    vec3 peak_c  = c + r * vec3(-1.2, 1.4, 0.0);
    vec3 eyeL_c  = c + r * vec3(-0.7, 1.6, -0.95);
    vec3 eyeR_c  = c + r * vec3(-0.7, 1.6, 0.95);
    vec3 hair1_c = c + r * vec3(-0.5, 2.6, 0.0);
    vec3 hair2_c = c + r * vec3(-0.45, 2.7, -0.2);
    vec3 hair3_c = c + r * vec3(-0.45, 2.7, 0.2);

    // Compute SDF for each part
    float head = sdfSphere(local_p, head_c, 0.9 * r);
    float body = sdfSphere(local_p, body_c, 1.3 * r);
    float tail = sdfEllipsoid(local_p, tail_c, vec3(1.2 * r, 0.6 * r, 0.6 * r));
    vec3 local_p_wing = rotateXYZ(local_p, wingL_c, vec3(0.0, 0.0, -20.0));
    float wingL = sdfEllipsoid(local_p_wing, wingL_c, vec3(0.6 * r, 1.0 * r, 1.0 * r));
    float wingR = sdfEllipsoid(local_p_wing, wingR_c, vec3(0.6 * r, 1.0 * r, 1.0 * r));
    float peak = sdfEllipsoid(local_p, peak_c, vec3(0.25 * r, 0.4 * r, 0.2 * r));
    float eyeL = sdfSphere(local_p, eyeL_c, 0.15 * r);
    float eyeR = sdfSphere(local_p, eyeR_c, 0.15 * r);

    vec3 local_p_hair = rotateXYZ(local_p, hair1_c, vec3(0.0, 0.0, 30.0));
    float hair1 = sdfEllipsoid(local_p_hair, hair1_c, vec3(0.08, 0.3, 0.2) * r * 1.7);
    local_p_hair = rotateXYZ(local_p, hair1_c, vec3(0.0, 45.0, 30.0));
    float hair2 = sdfEllipsoid(local_p_hair, hair2_c, vec3(0.08, 0.3, 0.2) * r * 1.7);
    local_p_hair = rotateXYZ(local_p, hair1_c, vec3(0.0, -45.0, 30.0));
    float hair3 = sdfEllipsoid(local_p_hair, hair3_c, vec3(0.08, 0.3, 0.2) * r * 1.7);

    // Combine the parts smoothly
    body = sdfUnionSmooth(sdfUnionSmooth(head, body, .1), tail, .1);
    float wings = sdfUnion(wingL, wingR);
    body = sdfUnion(body, peak);
    body = sdfSubtraction(body, eyeL);
    body = sdfSubtraction(body, eyeR);
    
    if (birb_hair) {
        float hair = sdfUnion(sdfUnion(hair1, hair2), hair3);
        body = sdfUnionSmooth(body, hair, .01);
    }   
    return sdfUnion(body, wings);
}

float sdfRiver(vec3 p)
{
    float riverbody = sdfBox(p, vec3(0.0, -0.65, 0.0), vec3(2.0, 0.1, 100.0));
    float t = time_cycle;
    float wave = 0.0015 * cos(8.0 * p.z - t * 15.0); // Curvy wave effect
    return riverbody + wave;
}

float sdfTree(in vec3 p, vec3 c)
{
    return sdfCone(p, c + vec3(0.0, 1.0, 0.0), 0.18, 2.0);
    // float leaves = sdfCone(p, c + vec3(0.0, 1.2, 0.0), 0.2, 0.3);
    // float leaves2 = sdfCone(p, c + vec3(0.0, 1.08, 0.0), 0.2, 0.4);
    // float leaves3 = sdfCone(p, c + vec3(0.0, 0.96, 0.0), 0.2, 0.5);
    // float trunk = sdfVerticleCapsule(p, c - vec3(0.0, 2.0, 0.0), 0.03, 3.0);
    // leaves = sdfUnion(leaves, leaves2);
    // leaves = sdfUnion(leaves, leaves3);
    // return sdfUnion(leaves, trunk);
}

float sdfGate(in vec3 p, vec3 c)
{
    p -= c;
    float r = 0.05;
    float pillarL = sdfVerticleCapsule(p, vec3(-0.7, -1.0, 0.0), r, 3.05);
    float pillarR = sdfVerticleCapsule(p, vec3(0.7, -1.0, 0.0), r, 3.05);
    // float pillarC = sdfVerticleCapsule(p, c + vec3(0.0, 1.7, 0.0), r, 0.3);
    float pillarC = sdfBox(p, vec3(0.0, 1.85, 0.0), vec3(0.08, 0.15, 0.05));
    float barTop = sdfBox(p, vec3(0.0, 2.0, 0.0), vec3(1.2, 0.06, 0.07));
    float barTop2 = sdfBox(p, vec3(0.0, 2.08, 0.0), vec3(1.3, 0.06, 0.10));
    float barBottom = sdfBox(p, vec3(0.0, 1.7, 0.0), vec3(1.05, 0.05, 0.07));

    float s = sdfUnion(pillarL, pillarR);
    s = sdfUnion(s, pillarC);
    s = sdfUnion(s, barTop);
    s = sdfUnion(s, barTop2);
    s = sdfUnion(s, barBottom);
    return s;
}

float sdfBribs(in vec3 p)
{
    float birb1 = 1000.0;
    float birb2 = 1000.0;
    float birb3 = 1000.0;
    if (p.z < 2.0) {
        vec3 birb_move = vec3(
            -0.3,
            0.03 * sin(iTime / 2.0),
            0.0
        );
        float birb_rotate = 60.0 * time_cycle;

        // White birb
        birb1 = sdfBorb(p, vec3(0.3, -0.6, -1.3 + 1.) + birb_move, -80.0, true);
        // Pink birb
        birb2 = sdfBorb(p, vec3(-0.5, -0.5, -0.3 + 1.) + birb_move, 160.0 + birb_rotate, true);
        // Yellow birb
        birb3 = sdfBorb(p, vec3(0.6, -0.5, 0.6 + 1.) + birb_move, 40.0 - birb_rotate * 2.0 , true);
    }
    return sdfUnion(birb1, sdfUnion(birb2, birb3));
}

float sdfFog(in vec3 p)
{
    p = vec3(p.x, p.y, p.z);

    // Bend the volume up based on p.x
    // p.y -= abs(p.x) * 0.1 - 2.0;

    float fog1 = sdfBox(p, vec3(0.0, -3.0, 10.0), vec3(500.0, 2.7, 500.0));
    float fog2 = sdfBox(p, vec3(0.0, -4.0, 20.0), vec3(500.0, 4.0, 8.0));
    float fog = sdfUnion(fog1, fog2);

    return fog;

    // float cloud1 = sdfSphere(p, vec3(0.0, 2.0, 5.0), 5.5);
    // float cloud2 = sdfSphere(p, vec3(0.0, 2.0, 5.0), 5.5);
    // float cloud = sdfUnion(cloud1, cloud2);
    // return sdfUnion(fog, cloud);
}


float sdfScene(vec3 p, bool record_hit)
{
    
    float background = sdfBox(p, vec3(0.0, 0.0, 50.0), vec3(100.0, 100.0, 1.0));
    // March background only if high enough
    if (p.y > 10.0) {
        hit_id.id = 5; // Background
        return background;
    }

    // Calculate the SDF for each 5 objects
    float ground = sdfCurvyGround(p, -0.1);
    float mountain1 = sdfOctahedron(p, vec3(-8.0, -2.5, 25.0), 4.0);
    float mountain2 = sdfOctahedron(p, vec3(-7.0, -2.5, 30.0), 5.0);
    float mountain3 = sdfOctahedron(p, vec3(-11.0, -2.5, 35.0), 7.0);
    float mountain4 = sdfOctahedron(p, vec3(-13.0, -2.5, 27.0), 7.5, vec3(10.0, 0.0, 0.0));  
    float mountain5 = sdfOctahedron(p, vec3(8.0, -2.5, 25.0), 4.0);
    float mountain6 = sdfOctahedron(p, vec3(7.0, -2.5, 30.0), 5.0);
    float mountain7 = sdfOctahedron(p, vec3(14.0, -2.5, 35.0), 7.0);
    float mountain8 = sdfOctahedron(p, vec3(8.0, -2.5, 20.0), 4.0);

    ground = sdfUnionSmooth(ground, mountain1, 0.5);
    ground = sdfUnionSmooth(ground, mountain2, 0.5);
    ground = sdfUnionSmooth(ground, mountain3, 0.5);
    ground = sdfUnionSmooth(ground, mountain4, 0.5);
    ground = sdfUnionSmooth(ground, mountain5, 0.5);
    ground = sdfUnionSmooth(ground, mountain6, 0.5);
    ground = sdfUnionSmooth(ground, mountain7, 0.5);
    ground = sdfUnionSmooth(ground, mountain8, 0.5);


    // change p for riverbed with respect to p.z
    vec3 riverbed_p = p + vec3(sin(p.z * 0.3), 0.0, 0.0);
    float riverbed = sdfBox(riverbed_p, vec3(0.0, 0.0, 0.0), vec3(1.5, 0.8, 100.0));
    ground = sdfSubtractionSmooth(ground, riverbed, 0.5);
    float riverbody = sdfRiver(p);


    // Distance Culling for birbs
    float birb1 = 1000.0;
    float birb2 = 1000.0;
    float birb3 = 1000.0;
    if (p.z < 2.0) {
        vec3 birb_move = vec3(
            -0.3,
            0.03 * sin(iTime / 2.0),
            0.0
        );
        float birb_rotate = 60.0 * time_cycle;

        // White birb
        birb1 = sdfBorb(p, vec3(0.3, -0.6, -1.3 + 1.) + birb_move, -80.0, true);
        // Pink birb
        birb2 = sdfBorb(p, vec3(-0.5, -0.5, -0.3 + 1.) + birb_move, 160.0 + birb_rotate, true);
        // Yellow birb
        birb3 = sdfBorb(p, vec3(0.6, -0.5, 0.6 + 1.) + birb_move, 40.0 - birb_rotate * 2.0 , true);
    }


    // float sun_move = float(int(iTime * 60.0) % 600) / 240.0;
    float sun_move = time_cycle * 4.0;
    float sun = sdfSphere(p, vec3(1.3, sun_move - 3.2, 50.0), 2.0);

    // Tree line
    float tree = 1000.0;
    // Distance Culling for trees
    if (p.z > 5.0 && p.z < 15.0) {
        vec3 tree_c = vec3(-13., -0.2, 10.0);
        // Generate positions for multiple trees
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 28; j++) {
                if (hash(float(i + j)) > 0.9) {
                    continue; // Randomly skip some trees
                }
                // xz direction displacement
                vec3 curr_c = tree_c + vec3(float(j) * 0.9 + float(i % 2) * 0.33, 
                                            0.0,
                                            float(i) * 2.5);
                // Avoid intersecting with the river
                if (curr_c.x > -1.2 && curr_c.x < 1.2) {
                    continue;
                }
                // y direction displacement
                curr_c.y += pow(curr_c.x, 2.0) * 0.007 + float(i) * 0.10;
                tree = sdfUnion(tree, sdfTree(p, curr_c));
            }
        }
    }

    // Gate
    float gate = 1000.0;
    if (p.z < 11.0 && p.x < 1.5 && p.x > -1.5) {
        gate = sdfGate(p, vec3(0.0, -0.4, 10.0));
    }

    // Combine the SDF for all objects
    float objects[] = float[](
        ground,
        birb1,
        birb2, 
        birb3,
        background,
        sun,
        riverbody,
        tree,
        gate
    );
    // Assign object ids for coloring
    int object_ids[] = int[](
        1,
        2,
        3, 
        4,
        5,
        6,
        7,
        8,
        9
    );
    float s = 1000.0; // set a large initial distance for union
    for (int i = 0; i < objects.length(); i++) {
        s = sdfUnion(s, objects[i]);
        // Record the closest object hit
        if (record_hit && s < hit_id.dist) {
            hit_id.dist = s;
            hit_id.id = object_ids[i];
        }
    }

    return s;
}

/** 
 * Overload sdfScene without hit_id update
 * E.g. we don't need to know what the object is in normal calculation
 */
float sdfScene(vec3 p)
{
    bool record_hit = true; // TODO: Should be false,
                            // but if I disable hit_id in normal calculation, 
                            // there will be artifacts in reflection... Don't know why yet
    return sdfScene(p, record_hit);
}

/**
 * Normal calculation for SDF Scene
 * @param p: intersection point query
 * @return norm: normal at the intersection point
 */
vec3 normal(vec3 p)
{
    float s = sdfScene(p); // sdf value in p
    float dx = 0.011; // Can't get smaller than this.. 
                      // 0.010 will produce wired artifacts in reflection... Don't know why yet

    vec3 norm = vec3(
        sdfScene(p + vec3(dx, 0.0, 0.0)) - s, // dsx
        sdfScene(p + vec3(0.0, dx, 0.0)) - s, // dsy
        sdfScene(p + vec3(0.0, 0.0, dx)) - s  // dsz
    );
    return normalize(norm);
}

float rayMarching(vec3 origin, vec3 dir)
{ 
    float s = 0.0; // distance
    for(int i = 0; i < 200; i++)
    {
        vec3 p = origin + dir * s;
        float dist = sdfScene(p, true); // sdf value in p
        s += dist; // update the distance
        if (s > 100.0 || abs(dist) < 0.002) {
            break;
        }
    }
    
    return s;
}

float rayMarchingBirbsOnly(vec3 origin, vec3 dir)
{ 
    float s = 0.0; // distance
    for(int i = 0; i < 128; i++)
    {
        vec3 p = origin + dir * s;
        float dist = sdfBribs(p);
        s += dist; // update the distance
        if (s > 20.0 || abs(dist) < 0.002) {
            break;
        }
    }
    
    return s;
}

vec3 sunrise(float t) {
    t = pow(t, 0.4);
    t = clamp(t, 0.0, 1.0);

    // Manual 8-stop gradient (could be put in a uniform array)
    const vec3 cols[8] = vec3[8](
        vec3(0.04, 0.10, 0.29),
        vec3(0.08, 0.18, 0.42),
        vec3(0.15, 0.29, 0.60),
        vec3(0.37, 0.35, 0.65),
        vec3(0.77, 0.43, 0.25),
        vec3(0.97, 0.62, 0.30),
        vec3(1.00, 0.83, 0.48),
        vec3(1.00, 0.96, 0.82)
    );

    // Scale t so each segment spans equal length (7 segments here)
    float f = 0.1 + t * 3.8;
    int i  = int(f);
    float w = fract(f);

    return mix(cols[i], cols[i+1], w);  // linear interpolate
}

vec3 stars(in vec3 p)
{
    vec3 c = vec3(0.);
    float res = iResolution.y*1.2;
    // float res = 1500.0;

	for (float i=0.;i<4.;i++)
    {
        vec3 q = fract(p*(.15*res))-0.5;
        vec3 id = floor(p*(.15*res));
        vec2 rn = hash33(id).xy;
        float c2 = 1.-smoothstep(0.,0.8,length(q));
        c2 *= step(rn.x,.0005+i*i*0.001);
        c += c2*(mix(vec3(1.0,0.49,0.1),vec3(0.75,0.9,1.),rn.y)*0.1+0.9);
        p *= 1.3;
    }
    return c*c*.8;
}

vec3 starTrail(in vec3 p)
{
    // transform p to rotate the stars
    float t = iTime * 2.0;
    vec3 color = vec3(0.0);
    
    for (int i = 0; i < 10; i++) {
        vec3 p = rotateXYZ(p, vec3(0.0), vec3(0.0, 0.0, t - float(i) * 0.8));
        color += stars(p) * pow(0.85, float(i));
    }

    return color;
}

// ----------------------------------------------------------------------------------
// Shading
// ----------------------------------------------------------------------------------
/**
 * A Copy of phongShading for handling reflection.
 * Because GLSL does not support recursive function calls.
 * This is with further reflection part removed because we just need one bounce.
 */
vec3 phongShadingReflection(
    vec3 p,                  // fragment position
    vec3 n,                  // surface normal
    vec3 ray_dir,            // viewing ray direction
    vec3 origin,             // camera/viewer position
    Light light,             // point light
    vec3 sunDir,             // additional sun light direction
    float sunDiffuse,        // sun diffuse strength
    float brightness_scale   // sunrise/sunset scaling
) {
    vec3 l = normalize(light.pos - p);
    vec3 eye = normalize(origin - p);
    vec3 r = reflect(-l, n);

    float amb = light.ambientStrength;
    float dif = max(dot(n, l), 0.0) * light.diffuseStrength;
    float spec = pow(max(dot(r, eye), 0.0), light.shininess) * light.specularStrength;
    float sunDif = max(dot(n, sunDir), 0.0) * sunDiffuse;

    vec3 color = vec3(1.0);
    float birb_brightness = 1.2;

    // Material colors by ID
    switch (hit_id.id) {
        case 1: // Ground
            color = vec3(0.8, 0.86, 0.88);
            if (p.y > 2.0) {
                color *= mix(1.0, 1.2, p.y - 2.0);
            }
            break;
        case 2: // Birb 1
            color = vec3(1.0, 0.89, 0.97) * birb_brightness;
            break;
        case 3: // Birb 2
            color = vec3(1.0, 0.61, 0.78) * birb_brightness;
            break;
        case 4: // Birb 3
            color = vec3(0.99, 0.79, 0.68) * birb_brightness;
            break;
        case 5: { // Background Sky
            vec3 color1 = sunrise(time_cycle + 0.2);
            vec3 color2 = sunrise(time_cycle);
            color = mix(color1, color2, (p.y + 20.0) / 50.0) * brightness_scale;

            // add some stars
            if (p.y > 0.0) {
                vec3 star = starTrail(ray_dir);
                color += star * (1.0 - time_cycle);
            }
            return color;
        }
        case 6: // Sun
            vec3 sun_color = vec3(1.0, 0.88, 0.57);
            return sun_color;
        
        case 7: { // River
            return vec3(0.79, 0.89, 1.0);
        }

        case 8: // Tree
            color = vec3(0.65, 0.77, 0.54);
            float disp = pow(p.x, 2.0) * 0.007 + 0.15;
            if (p.y > disp) {
                color = mix(color, vec3(1.0) * 1.5, p.y - disp);
            }
            break;

            
        case 9: // Gate
            color = vec3(0.76, 0.55, 0.55);
            if (p.y > 1.29 && p.y < 1.35) {
                color = mix(color, vec3(1.0), (p.y - 1.20) * 10.0);
            }
            if (p.y > 1.55) {
                color = mix(color, vec3(1.0), (p.y - 1.55) * 10.0);
            } 
            break;

        default:
            return vec3(0.13, 1.0, 0.0); // unexpected, debug color
    }

    // Shadow computation, ignore in reflection to save resources
    // float s = rayMarching(p + n * 0.02, l);
    // if (s < length(light.pos - p)) {
    //     dif *= 0.2;
    // }

    // Final shaded color
    return (amb + dif + spec + sunDif) * brightness_scale * color * light.color;

}

vec3 phongShading(
    vec3 p,                  // fragment position
    vec3 n,                  // surface normal
    vec3 ray_dir,            // viewing ray direction
    vec3 origin,             // camera/viewer position
    Light light,             // point light
    vec3 sunDir,             // additional sun light direction
    float sunDiffuse,        // sun diffuse strength
    float brightness_scale   // sunrise/sunset scaling
) {
    vec3 l = normalize(light.pos - p);
    vec3 eye = normalize(origin - p);
    vec3 r = reflect(-l, n);

    float amb = light.ambientStrength;
    float dif = max(dot(n, l), 0.0) * light.diffuseStrength;
    float spec = pow(max(dot(r, eye), 0.0), light.shininess) * light.specularStrength;
    float sunDif = max(dot(n, sunDir), 0.0) * sunDiffuse;

    vec3 color = vec3(1.0);
    float birb_brightness = 1.2;

    // Material colors by ID
    switch (hit_id.id) {
        case 1: // Ground
            color = vec3(0.8, 0.86, 0.88);
            if (p.y > 2.0) {
                // color *= mix(1.0, 1.2, p.y - 2.0);
                color = mix(color, vec3(1.0), p.y - 2.0);
            }
            break;

        case 2: // Birb 1
            color = vec3(1.0, 0.89, 0.97) * birb_brightness;
            break;

        case 3: // Birb 2
            color = vec3(1.0, 0.61, 0.78) * birb_brightness;
            break;

        case 4: // Birb 3
            color = vec3(0.99, 0.79, 0.68) * birb_brightness;
            break;

        case 5: // Background Sky
            vec3 color1 = sunrise(time_cycle + 0.2);
            vec3 color2 = sunrise(time_cycle);
            float a = smoothstep(-1000.0, 2000.0, pow(p.y, 2.0) * 4.0 + pow(p.x, 2.0));
            color = mix(color1, color2, a) * brightness_scale;
            // color = mix(color1, color2, (p.y + 20.0) / 50.0) * brightness_scale;

            // add some stars
            if (p.y > 0.0) {
                vec3 star = starTrail(ray_dir);
                color += star * (1.0 - time_cycle) * 1.25;
            }
            return color;

        case 6: // Sun
            vec3 sun_color = vec3(1.0, 0.88, 0.57) * 1.1;
            return sun_color;

        case 7: // River
            vec3 water_color = vec3(0.79, 0.89, 1.0);

            vec3 reflect_dir = reflect(ray_dir, n);
            float reflect_s = rayMarching(p + n * 0.01, reflect_dir);
            vec3 reflect_p = p + reflect_dir * reflect_s;
            vec3 reflect_n = normal(reflect_p);
            vec3 reflect_color = phongShadingReflection(
                reflect_p,
                reflect_n,
                reflect_dir,
                CAM_POS,
                light,
                sunDir,
                sunDiffuse,
                brightness_scale
            );

            color = mix(water_color, reflect_color, 0.7);
            break;

        case 8: // Tree
            color = vec3(0.65, 0.77, 0.54);
            float disp = pow(p.x, 2.0) * 0.007 + 0.15;
            if (p.y > disp) {
                color = mix(color, vec3(1.0) * 1.5, p.y - disp);
            }
            break;

        case 9: // Gate
            color = vec3(0.76, 0.55, 0.55);
            if (p.y > 1.29 && p.y < 1.35) {
                color = mix(color, vec3(1.0), (p.y - 1.20) * 10.0);
            }
            if (p.y > 1.55) {
                color = mix(color, vec3(1.0), (p.y - 1.55) * 10.0);
            } 
            break;

        default:
            return vec3(0.13, 1.0, 0.0); // unexpected, debug color
    }

    // Shadow computation
    float s = rayMarching(p + n * 0.02, l);
    if (s < length(light.pos - p)) {
        dif *= 0.2;
    }

    // Fog effect
    // float fog = 1.0 - exp(-0.03 * p.z);
    // vec3 fog_color = vec3(0.49, 0.44, 0.4);

    // Final shaded color
    vec3 lit_color = (amb + dif + spec + sunDif) * brightness_scale * color * light.color;
    return  lit_color; // mix(lit_color, fog_color, fog);
}

// ----------------------------------------------------------------------------------
// Volumetric Fog Rendering
// ----------------------------------------------------------------------------------


float fogDensity(in vec3 p) 
{
    // p = p + vec3(time_cycle * 1.0, 0.0, time_cycle * 5.0);
    p.x += mod(iTime * 0.5, 100.0);
    vec3 q = (p - vec3(0.0,5.0,10.0)) * 3.0 - vec3(0.0,0.5,1.0)*iTime * 0.2;
    float f = fbm(q);
    float tor = 1.0 - sdfFog(p * 4.0) + f * 2.0;
    return clamp(tor * 0.08, 0.0, 1.0);
}

#define SHADOW_STEPS 1
#define SHADOW_LENGTH 1.0

float calcShadow(vec3 pos, Light light)
{

    vec3 lightDir = normalize(light.pos - pos);
    float shadow = 0.0;
    float step = SHADOW_LENGTH / float(SHADOW_STEPS);
    for (int i = 0; i < SHADOW_STEPS; i++) {
        pos += lightDir * step;
        shadow += fogDensity(pos);
    }
    return shadow;
}

vec3 fogShading(
    in vec3  pos,
    in Light light,
    in float density,
    in float alpha,
    out float newAlpha
) {
    // vec3 lightDir = normalize(light.pos - pos);
    // 1) shadow attenuation
    float shadow = calcShadow(pos, light);
    float atten  = exp(-shadow / float(SHADOW_STEPS) * 3.0);
    // 2) scattering
    vec3 color = atten * density * alpha * light.color;

    // 3) sky glow
    vec3 color_sky = vec3(1.0);
    // vec3 color_sky = sunrise(time_cycle);
    float skyT = exp(-fogDensity(pos) * 0.8);
    color += skyT * density * alpha * color_sky * light.color;

    // Check if the point is in shadow
    if (pos.z < 4.0) {
        HitID hit_id_copy = hit_id;
        float lightDist = length(light.pos - pos);
        float sceneDist = rayMarchingBirbsOnly(pos, normalize(light.pos - pos));
        if (sceneDist < lightDist - 0.001) {
            color *= 0.6; // darken the color
        } else {
            color *= 1.0; // brighten the color
        }
        hit_id = hit_id_copy;
    }

    // 4) update alpha
    newAlpha = alpha * (1.0 - density);
    return color;
}

#define MAX_STEPS     48 // 12
#define VOLUME_LENGTH 20.0 // 6.0

vec4 rayMarchVolume(vec3 origin, vec3 ray, Light light, float brightness_scale) 
{
    vec4 sum = vec4(0,0,0,1);
    float stepLen = VOLUME_LENGTH / float(MAX_STEPS);
    vec3 pos = origin + ray * jitter * stepLen;
    // vec3 pos = origin + ray * stepLen;

    // March ray into the scene
    float sceneHitDist = rayMarching(origin, ray);

    float contribScale = 1.0;
    for (int i = 0; i < MAX_STEPS; i++) {
        if (sum.a < 0.1) break;

        float d = clamp(fogDensity(pos) * contribScale, 0.0, 1.0);
        if (d > 0.001) {
            float newAlpha;
            vec3  contrib = fogShading(pos, light, d, sum.a, newAlpha);
            sum.rgb += contrib;
            sum.a    = newAlpha;
        }
        pos += ray * stepLen;

        // if hit the scene, break
        if (sceneHitDist < length(pos - origin)) {
            break;
        }

        // use larger step length to speed up
        if (i >= MAX_STEPS / 2) {
            if (contribScale <= 2.0) {
                contribScale *= 2.0;
                stepLen *= 2.0;
            }
        }
    }

    // Scene hit
    vec3 sceneHitPos = origin + ray * sceneHitDist;
    vec3 norm = normal(sceneHitPos);
    vec3 sunDir = vec3(0.0, 1.0, -1.0);
    // if (hit_id.id == 9) {
    //     sunDir = vec3(0.0, 1.0, 1.0);
    // }
    vec3 color = phongShading(
        sceneHitPos, norm, ray, origin,
        light, sunDir, 0.2,
        brightness_scale
    );
    sum.rgb = color * sum.a + sum.rgb * min(brightness_scale * 0.5, 0.9);
    sum.a = 1.0; // set alpha to 1.0 to make the scene opaque

    // }
    return sum;
}

vec4 postProcessing(vec2 uv, in vec4 color)
{   
    // Apply a vignette effect
    float strength = 0.6;
    float vignette = smoothstep(0.0, 1.0, length(uv) * 1.0);
    color.rgb = mix(color.rgb, vec3(0.0), vignette * strength);

    // vec3 col = color.rgb;
	// col *= 0.5 + 0.5*pow( 16.0*uv.x*uv.y*(1.0-uv.x)*(1.0-uv.y), 0.1 );

    // col = clamp(col,0.0,1.0);
    // col = col*0.6 + 0.4*col*col*(3.0-2.0*col) + vec3(0.0,0.0,0.04);
    // color.rgb = col;

    // Increase contrast
    float contrast = 1.3;
    color.rgb = (color.rgb - 0.5) * contrast + 0.5;

    // Brightness adjustment
    float brightness = 0.1;
    color.rgb += brightness;

    return color;
}


// On-Screen snow effect https://www.shadertoy.com/view/ldsGDn
#define SNOW_LAYERS 66

#define SNOW_DEPTH1 .3
#define SNOW_WIDTH1 .4
#define SNOW_SPEED1 .6

#define SNOW_DEPTH2 .1
#define SNOW_WIDTH2 .3
#define SNOW_SPEED2 .1

float snowing(in vec2 uv)
{
    const mat3 p = mat3(13.323122, 23.5112, 21.71123,
                        21.121200, 28.7312, 11.9312,
                        21.811200, 14.7212, 61.3934);
    vec2 mp = vec2(0.0);
    uv.x += mp.x * 4.0;    
    mp.y *= 0.25;
    float depth = smoothstep(SNOW_DEPTH1, SNOW_DEPTH2, mp.y);
    float width = smoothstep(SNOW_WIDTH1, SNOW_WIDTH2, mp.y);
    float speed = smoothstep(SNOW_SPEED1, SNOW_SPEED2, mp.y);
    float acc = 0.0;
    float dof = 5.0 * sin(iTime * 0.1);
    for (int i=0; i < SNOW_LAYERS; i++)
    {
        float fi = float(i);
        vec2 q = uv * (1.0 + fi * depth);
        float w = width * mod(fi * 7.238917, 1.0) - width * 0.1 * sin(iTime * 2. + fi);
        q += vec2(q.y * w, speed * iTime / (1.0 + fi * depth * 0.03));
        vec3 n = vec3(floor(q), 31.189 + fi);
        vec3 m = floor(n) * 0.00001 + fract(n);
        vec3 mp = (31415.9 + m) / fract(p * m);
        vec3 r = fract(mp);
        vec2 s = abs(mod(q, 1.0) - 0.5 + 0.9 * r.xy - 0.45);
        s += 0.01 * abs(2.0 * fract(10. * q.yx) - 1.); 
        float d = 0.6 * max(s.x - s.y, s.x + s.y) + max(s.x, s.y) - .01;
        float edge = 0.05 + 0.05 * min(.5 * abs(fi - 5. - dof), 1.);
        acc += smoothstep(edge, -edge, d) * (r.x / (1. + .02 * fi * depth));
    }
    return acc;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    jitter = hash(uv.x + uv.y * 57.0 + iTime * 0.1);

    time_cycle = getSyncedTimeCycle();

    float cam_move = time_cycle;

    vec3 origin = CAM_POS 
                + vec3(-0.2, -0.8, 1.7 - 1.0) 
                + vec3(cam_move * 1.0, cam_move * 0.0, -cam_move * 0.0);

    vec3 dir = normalize(vec3(uv.x, uv.y, 1.0));

    float brightness_scale = 1.4 + time_cycle * 0.3;

    Light light;
    light.pos = vec3(0.0, smoothstep(1.0, 20.0, time_cycle), 32.0);
    light.color = mix(sunrise(time_cycle), vec3(1.0), 0.2);
    light.ambientStrength = 0.5;
    light.diffuseStrength = 0.6;
    light.specularStrength = 0.7;
    light.shininess = 64.0;

    vec4 color = rayMarchVolume(origin, dir, light, brightness_scale);

    vec2 uv_snow = uv * iResolution.x / 1000.0 + vec2(cam_move * 0.5, 0.0);
    float snow = snowing(uv_snow) * brightness_scale * 0.4;
    
    color += vec4(vec3(snow) * mix(sunrise(time_cycle), vec3(1.0), 0.5), 0.0);

    color = postProcessing(uv, color);

    fragColor = vec4(color.rgb, 1.0); // Set the alpha to 1.0 for the final color
}

void main() 
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}