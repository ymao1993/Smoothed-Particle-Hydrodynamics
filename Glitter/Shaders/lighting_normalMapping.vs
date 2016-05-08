#version 430 core

//vertex attributes
layout(location = 0) in vec3 position;

//uniforms
layout(location = 0) uniform float radius;
layout(location = 1) uniform mat4 m2w;
layout(location = 2) uniform mat4 w2v;
layout(location = 3) uniform mat4 persp_proj;

out vec3 vpos; //view space position

void main()
{
   //set point size according to eye-point distance
   vec4 pos = w2v * m2w * vec4(position,1);
   vpos = pos.xyz / pos.w;
   float dist = length(vpos);
   gl_PointSize = radius / dist;
   gl_Position = persp_proj * pos;
}

