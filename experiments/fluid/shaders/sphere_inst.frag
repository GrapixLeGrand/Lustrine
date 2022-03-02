#version 430

//a tutorial https://paroj.github.io/gltut/Illumination/Tutorial%2013.html
out vec4 out_color;

uniform mat4 p;
uniform vec3 light_direction;

in vec3 v2f_normal;
in vec2 v2f_uv;
in vec3 v2f_position_view;
in vec3 v2f_position;
in vec4 selected_color;

vec3 ambiant = vec3(0.8, 0.1, 0.2);

float radius = 0.5;

void main()
{   

    vec2 uv = 2.0 * (v2f_uv - vec2(0.5, 0.5));
    float len = dot(uv, uv);

    if (len <= 1.0) {
        
        //uv.y *= -1;
        vec3 normal = vec3(uv, sqrt(1.0 - len)); //in camera position
        vec3 position = normalize(normal) * radius + v2f_position_view;
        vec4 clipPos = p * vec4(position, 1.0);
        float depth_ndc = clipPos.z / clipPos.w; 
        gl_FragDepth = ((gl_DepthRange.diff * depth_ndc) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;
    
        vec3 light_dir = light_direction;

        float diff = max(dot(normal, light_dir), 0.0);
        //vec3 reflect_dir = reflect(-light_dir, normal);
        //float spec = pow(max(dot(-position, reflect_dir), 0.0), shininess);
        out_color = selected_color * vec4(vec3(0.5), 1.0) + diff * selected_color;
    
    } else {
        discard;
    }

    
}