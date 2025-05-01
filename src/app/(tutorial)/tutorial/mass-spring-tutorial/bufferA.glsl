precision highp float;


#define COLOR
#define TUG // comment out to snap to cursor

varying vec2 vUv; // UV (screen) coordinates in [0,1]^2

uniform float iTime;
uniform float iTimeDelta;
uniform float iFrame;
uniform vec2 iResolution;
uniform vec4 iMouse;
uniform sampler2D iChannel0;


const float k = 0.5;
const float b = 0.0;
const float l = 0.5;
const float m = 1.0;

const vec2 g = vec2(0, -9.8) * 2.0;

const float[] bs = float[](b, b, b);
const float[] ks = float[](k, k, k);
const float[] ls = float[](l, l, l);
const float[] ms = float[](m, m, m);

const vec4[] bodies = vec4[](vec4(1,0.5,-7,15), vec4(0.3,1.5,-3,0), vec4(-0.5,2,15,4));
const int n = 3;

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

// Simulation

void main()
{
    vec2 fragCoord = vUv * iResolution.xy; // Calculate pixel coordinates from UV and resolution
    if (fragCoord.y < float(n) && fragCoord.x < 2.0) 
    {
    	if (iFrame < 1.0 || resolutionChanged(iChannel0, iResolution))
    	{
    	    // Initial conditions.
    	    gl_FragColor = bodies[int(fragCoord.y)];
    	}
        else if (fragCoord.x < 1.0)
        {
            int i = int(fragCoord.y);
            float t = min(iTimeDelta, 0.2);
            
            vec4 body = getBody(i, iChannel0);
            vec4 anchor = getBody(i-1, iChannel0);
            
            if (i > 0)
            {
                anchor.zw += g*t;
            }
            body.zw += g*t;
            
            // Spring equation
            vec2 dx = anchor.xy - body.xy;
            float l = ls[i];
            float x = l - length(dx);
            
            vec2 d = normalize(dx);
            vec2 dv = anchor.zw - body.zw;
            float v = dot(dv, d);
            
            float k = ks[i];
            float b = bs[i];
            float m = ms[i];
            
            bool isMouseDown = iMouse.z > 0.;
            vec2 f = (-k*x + b*v) * d;
            body.zw += f / m;
            
            const float drag = 0.98;
            body.zw *= drag;

            
            if (i < n-1)
            {
                k = ks[i+1];
                b = bs[i+1];
                l = ls[i+1];
                
                vec4 anchor = getBody(i+1, iChannel0);
                anchor.zw += g*t;
             
                // Spring equation
            	vec2 dx = body.xy - anchor.xy;
            	float x = l - length(dx);
            
            	vec2 d = normalize(dx);
            	vec2 dv = body.zw - anchor.zw;
            	//float v = dot(dv, d);
            
            	vec2 f = (k*x - b*v) * d;
            	body.zw += f / m;
            }
            
            // Integration
            body.xy += body.zw * t;
            
            if (isMouseDown && fragCoord.y > 2.0) {
                float scale = iResolution.y / 20.0;
    			vec2 middle = iResolution.xy / 2.0 + vec2(0, 5.0) * scale;
                
                vec2 mousePos = iMouse.xy;
                vec2 bodyPos = (body.xy * scale + middle);
                vec2 delta = mousePos - bodyPos;
            	
#ifdef TUG
                // tug
                body.zw += delta * 0.05;
#else
                // snap
                body.xy = (mousePos - middle) / scale;
                body.zw = vec2(0.,0.);
#endif
            }
            gl_FragColor = body;
        }
        else
        {
            int i = int(fragCoord.y);
            gl_FragColor = getBody(i, iChannel0);
        }
    }
    else
    {
        gl_FragColor = vec4(iResolution.x, iResolution.y,0., 0.);
    }
}
