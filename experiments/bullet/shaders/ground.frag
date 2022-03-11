#version 430

out vec4 out_color;

in vec2 v2f_tex_coord;
uniform sampler2D tex;

void main()
{
    vec2 c = v2f_tex_coord;
    out_color = texture(tex, 30.0 * c); //vec4(v2f_tex_coord, 0.0, 1.0);
}