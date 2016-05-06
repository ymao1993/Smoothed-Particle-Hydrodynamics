/*Modified from XRealityRender's Shader Utility Class*/

#ifndef SHADERUTILS_H
#define SHADERUTILS_H
#include <glad/glad.h>


/**
 * XRShaderUtils
 * Some utility functions for shader operations.
 *
 * @Author Yu Mao
 */
namespace ShaderUtils
{
		GLuint loadShader(const char * filename,
			GLenum shader_type = GL_FRAGMENT_SHADER,
			bool check_errors = true);

		GLuint loadShaderFromSrc(const char * source,
			GLenum shader_type,
			bool check_errors = true);

		GLuint linkShaderProgram(const GLuint * shaders,
			int shader_count,
			bool delete_shaders,
			bool binary = false,
			bool check_errors = true);

};

#endif