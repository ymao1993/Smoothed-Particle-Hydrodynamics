/*
* Copyright � 2012-2015 Graham Sellers
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice (including the next
* paragraph) shall be included in all copies or substantial portions of the
* Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*/

#include <glad/glad.h>
#include <cstdio>


namespace ShaderUtils
{
	GLuint loadShader(const char * filename, GLenum shader_type, bool check_errors)
	{
		GLuint result = 0;
		FILE * fp;
		size_t filesize;
		char * data;

		fp = fopen(filename, "rb");

		if (!fp) {
 		   printf("Cannot open file: %s\n", filename);
                   return -1;
		}

		fseek(fp, 0, SEEK_END);
		filesize = ftell(fp);
		fseek(fp, 0, SEEK_SET);

		data = new char[filesize + 1];

		if (!data)
			goto fail_data_alloc;

		fread(data, 1, filesize, fp);
		data[filesize] = 0;
		fclose(fp);

		result = glCreateShader(shader_type);

		if (!result)
			goto fail_shader_alloc;

		glShaderSource(result, 1, (const char**)&data, NULL);

		delete[] data;

		glCompileShader(result);

		if (check_errors)
		{
			GLint status = 0;
			glGetShaderiv(result, GL_COMPILE_STATUS, &status);

			if (!status)
			{
				char buffer[4096];
				glGetShaderInfoLog(result, 4096, NULL, buffer);
				printf("%s: %s\n", filename, buffer);
				goto fail_compile_shader;
			}
		}

		return result;

	fail_compile_shader:
		glDeleteShader(result);

	fail_shader_alloc:;
	fail_data_alloc:
		return result;
	}

	GLuint loadShaderFromSrc(const char * source,
		GLenum shader_type,
		bool check_errors)
	{
		GLuint sh;

		sh = glCreateShader(shader_type);

		const char * strings[] = { source };
		glShaderSource(sh, 1, strings, nullptr);

		glCompileShader(sh);

		if (check_errors)
		{
			GLint status = 0;
			glGetShaderiv(sh, GL_COMPILE_STATUS, &status);

			if (!status)
			{
				char buffer[4096];
				glGetShaderInfoLog(sh, 4096, NULL, buffer);
				printf("%s\n", buffer);
				goto fail_compile_shader;
			}
		}

		return sh;

	fail_compile_shader:
		glDeleteShader(sh);

		return 0;
	}


	GLuint linkShaderProgram(const GLuint * shaders,
		int shader_count,
		bool delete_shaders,
		bool binary,
		bool check_errors)
	{
		int i;

		GLuint program;

		program = glCreateProgram();

		for (i = 0; i < shader_count; i++)
		{
			glAttachShader(program, shaders[i]);
		}

		if (binary) 
			glProgramParameteri(program, GL_PROGRAM_BINARY_RETRIEVABLE_HINT, GL_TRUE);

		glLinkProgram(program);

		if (check_errors)
		{
			GLint status;
			glGetProgramiv(program, GL_LINK_STATUS, &status);

			if (!status)
			{
				char buffer[4096];
				glGetProgramInfoLog(program, 4096, NULL, buffer);
				glDeleteProgram(program);
				return 0;
			}
		}

		if (delete_shaders)
		{
			for (i = 0; i < shader_count; i++)
			{
				glDeleteShader(shaders[i]);
			}
		}

		return program;
	}
}
