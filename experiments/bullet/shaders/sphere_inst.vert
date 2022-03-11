#version 430

layout (location = 0) in vec3 vertex_position;
layout (location = 1) in vec2 vertex_uv;
layout (location = 2) in vec3 vertex_normal;
layout (location = 3) in vec3 world_position;
layout (location = 4) in vec4 color;

uniform mat4 vp; //only view projection, m is the I matrix
uniform mat4 view;
uniform mat3 view_inv;
uniform float scale;

out vec3 v2f_normal;
out vec2 v2f_uv;
out vec3 v2f_position_view;
out vec3 v2f_position;
out vec4 selected_color;


void main()
{   
    selected_color = color;
    
    v2f_uv = vertex_uv;
    v2f_normal = vertex_normal;
    vec4 temp = view * vec4(world_position, 1.0);
    v2f_position_view = temp.xyz / temp.w;
    v2f_position = world_position;

    vec3 pos = (view_inv * (vertex_position * scale) + world_position); //world position of the rotated vertex of the billboard

    gl_Position = vp * vec4(pos, 1.0);
}