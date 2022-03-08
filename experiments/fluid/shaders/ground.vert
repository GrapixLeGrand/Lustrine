#version 430

layout (location = 0) in vec3 vertex_position;
layout (location = 1) in vec2 vertex_uv;
layout (location = 2) in vec3 vertex_normal;

uniform mat4 mvp;

out vec2 v2f_tex_coord;

void main()
{
    v2f_tex_coord = vertex_uv;
    vec3 position = vertex_position;
    position.xz += 15.0;
    gl_Position = mvp * vec4(position, 1.0); 
}