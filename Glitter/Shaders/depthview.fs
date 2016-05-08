#version 430 core

//layout(binding = 0, r32f) readonly uniform image2D depth;
layout(location = 0) uniform sampler2D depth;

out vec4 color;

//
//transform depth from "squezzed" unit cube to normalized view space 
//
float LinearizeDepth(float depth)
{
    float zNear = 0.1;
    float zFar  = 100.0; 
    return (2.0 * zNear) / (zFar + zNear - depth * (zFar - zNear));
}

void main()
{
	float depth = texelFetch(depth, ivec2(gl_FragCoord.xy),0).x;
	//imageLoad(depth, ivec2(gl_FragCoord.xy)).x;

	depth = LinearizeDepth(depth);

	color = vec4(depth,depth,depth,1);
}