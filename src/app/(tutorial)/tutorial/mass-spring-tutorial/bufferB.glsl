#define COLOR

varying vec2 vUv; // UV (screen) coordinates in [0,1]^2


uniform float iTime;
uniform float iTimeDelta;
uniform float iFrame;
uniform vec2 iResolution;
uniform vec4 iMouse;

uniform sampler2D iChannel0;
uniform sampler2D iChannel1;

const vec4[] colors = vec4[](vec4(1,0,0,1), vec4(0,1,0,1), vec4(0,0,1,1));
const int n = 3;

float segmentDistance( vec2 p, vec2 a, vec2 b )
{
    vec2 pa = p - a, ba = b - a;
	float h = clamp(dot(pa,ba) / dot(ba,ba), 0., 1.);	
	return length(pa - ba * h);
}

vec3 toLinear(vec3 v) {
  	return pow(v, vec3(2.2));
}

vec3 tosRGB( vec3 color )
{
    return pow(color, vec3(1.0/2.2));
}

vec4 outputColor = vec4(vec3(0), 1);
vec4 color = vec4(1.0);

void drawLine( vec2 coords, vec2 p0, vec2 p1, float thickness )
{
    float sD = segmentDistance(coords, p0, p1);
	float a = 1.0 - clamp(sD - thickness / 2.0 + .5, 0.0, 1.0);
    
    outputColor = mix(outputColor, color, a * color.a);
}

void drawSpring( vec2 coords, vec2 p0, vec2 p1, float thickness, int loops )
{
    vec2 d = p1 - p0;
    if (length(d) < 0.001)
        return;
    
    vec2 dir = normalize(d);
    vec2 per = vec2(dir.y, -dir.x);
    
    vec2 st = d / float(loops * 2);
    vec2 last = p0 + per * thickness / 2.0 + st / 2.0;
    vec2 sw = -thickness * per;
    float th = 1.0;
    
    drawLine(coords, p0, last, th);
    
    for (int i=0; i<loops*2-1; i++)
    {
        vec2 next = last + st + sw;
        sw = -sw;
        drawLine(coords, last, next, th);
        last = next;
    }
    
    drawLine(coords, last, p1, th);
}

void drawCircle( vec2 coords, vec2 center, float radius )
{
    float sD = distance(coords, center);
	float a = 1.0 - clamp(sD - radius / 2.0 + .5, 0.0, 1.0);
    
        outputColor = mix(outputColor, color, a * color.a);
}

vec4 getBody( int i, sampler2D tex_sampler )
{
    if (i < 0)
    {
        return vec4(0.0); 
    }
        
    return texelFetch(tex_sampler, ivec2(0, i), 0);
}

vec4 getBodyPrevious( int i, sampler2D tex_sampler )
{
    if (i < 0)
    {
        return vec4(0.0); 
    }
        
    return texelFetch(tex_sampler, ivec2(1, i), 0);
}

bool resolutionChanged( sampler2D tex_sampler, vec2 iResolution )
{
    return iResolution != texelFetch(tex_sampler, ivec2(5, 5), 0).rg;
}


void main()
{
    vec2 fragCoord = vUv * iResolution.xy; // Calculate pixel coordinates from UV and resolution
    float scale = iResolution.y / 20.0;
    vec2 middle = iResolution.xy / 2.0 + vec2(0, 5.0) * scale;
    
    vec2 uv = fragCoord / iResolution.xy;
    
    outputColor = (resolutionChanged(iChannel1, iResolution) ? 0.0 : .998) * texture(iChannel1, uv);
    for (int i=0; i<n; i++)
    {
#ifdef COLOR
        color = colors[i];
#endif        
        drawLine(fragCoord, 
                 getBodyPrevious(i, iChannel0).xy * scale + middle, 
                 getBody(i, iChannel0).xy * scale + middle,
                 0.5);
    }
    
    gl_FragColor = outputColor;
    
    if (fragCoord.x < 6.0) 
    {
        gl_FragColor = vec4(iResolution, 0, 0);
    }
}