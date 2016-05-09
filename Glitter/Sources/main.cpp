// Local Headers
#include "glitter.hpp"

// System Headers
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// Standard Headers
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <math.h>
#include <algorithm>

#include "SOIL2.h"

#include "ShaderUtils.hpp"
#include "SPHSimulator.hpp"

#include "timer.hpp"

//application
GLFWwindow* mWindow;

//mouse input
static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

static void onRelease();
static void onPress();
static void duringPress();
static bool is_mouse_pressed = false;
static double cursor_pre_x;
static double cursor_pre_y;
static double cursor_cur_x;
static double cursor_cur_y;
static void arcball_update();
static const double CURSOR_POS_INVALID = -2;
static const double speed_factor = 2.;

//transformation matrices
static mat4 m2w;
static mat4 w2v;
static mat4 pers_proj;

//SPH Simulation
static SPHSim::BoxDef box;
static float delta = 1./60;
static SPHSim::SPHSimulator* simulator;

//Lighting
static vec3 ambient(0.2,0.0,0.0);
static vec3 diffuse(0.5,0.0,0.0);
static vec3 specular(3.0,3.0,3.0);
static vec3 lightDir(-0.5774,-0.5774,-0.5774);
static float shineness = 2;

//**SHADERS**//
#define PARTICLE_SPLAT_SIZE 400
//
//pointSplat
//
static const char* pointSplat_files[] = 
{
 "../Glitter/Shaders/pointSplat.vs",
 "../Glitter/Shaders/pointSplat.fs"
};
enum
{
  POINTSPLAT_LOCATION_POS,
};
enum 
{
   POINTSPLAT_LOCATION_R,
   POINTSPLAT_LOCATION_M2W,
   POINTSPLAT_LOCATION_W2V,
   POINTSPLAT_LOCATION_PERSP_PROJ
};
static GLuint pointSplat_program;
static GLuint pointSplat_vao;
static GLuint vbo_pointSplat_pos;

//
//gaussian blur
//
static const char* gaussianBlur_files[] = 
{
  "../Glitter/Shaders/GaussianBlur.cs"
};
static GLuint gaussianBlur_program;

//
//texture mapping
//
static const char* lighting_files[] = 
{
  "../Glitter/Shaders/lighting_normalMapping.vs",
  "../Glitter/Shaders/lighting_normalMapping.fs"
};
enum
{
   LIGHTING_LOCATION_POS,
};
enum 
{
   LIGHTING_LOCATION_R,
   LIGHTING_LOCATION_M2W,
   LIGHTING_LOCATION_W2V,
   LIGHTING_LOCATION_PERSP_PROJ,
   LIGHTING_LOCATION_AMBIENT,
   LIGHTING_LOCATION_DIFFUSE,
   LIGHTING_LOCATION_SPECULAR,
   LIGHTING_LOCATION_LDIR,
   LIGHTING_LOCATION_SHINESNESS,
};

static GLuint lighting_program;
static GLuint lighting_vao;
static GLuint vbo_lighting_pos;

//
//depth visualization
//
static const char* depthview_files[] =
{
 "../Glitter/Shaders/depthview.vs",
 "../Glitter/Shaders/depthview.fs",  
};
enum
{
  DEPTHVIEW_LOCATION_POS,
};
static GLuint depthview_program;
static GLuint depthview_vao;
static GLuint vbo_depthview_pos;

//
static GLuint fbo;
static GLuint normalBuffer;
static GLuint depthBuffer;
static GLuint blurredNormal;


//forward declarations
static void init();
static void render();
static void update();
static void onExit();

static void initShaders();
static void init_pointSplat();
static void init_map();
static void init_gaussianBlur();
static void init_depthview();
static void initSPH();
static void initCamera();
static void initOfflineRendering();
static void initPostprocessing();
static void renderBox();
static void renderParticles(bool offscreen);
static void renderParticlesAsPoints();
static void renderFluid();
static void blurTexture(GLuint input, GLuint output, int width, int height);
static void setTransform();
static void outputTexture2File(GLuint texture, const char* file, GLenum format, GLenum type, int channel);

//tests
static void testLoadImage();
static void visualizeDepth(GLuint depth);

int main(int argc, char * argv[]) {

    // Load GLFW and Create a Window
    glfwInit();
    mWindow = glfwCreateWindow(mWidth, mHeight, "Parallel Fluid Simulator", nullptr, nullptr);
 
    // Check for Valid Context
    if (mWindow == nullptr) {
        fprintf(stderr, "Failed to Create OpenGL Context");
        glfwTerminate();
        return EXIT_FAILURE;
    }

    //set mouse callback
    glfwSetMouseButtonCallback(mWindow, mouse_button_callback);
    glfwSetScrollCallback(mWindow, scroll_callback);

    // Create Context and Load OpenGL Functions
    glfwMakeContextCurrent(mWindow);
    gladLoadGL();
    fprintf(stderr, "OpenGL %s\n", glGetString(GL_VERSION));

    // Enable Depth Test
    glEnable(GL_DEPTH_TEST);

    // Initialization
    init();

    // Rendering Loop
    while (glfwWindowShouldClose(mWindow) == false) {
        
        //exit when press esc
        if (glfwGetKey(mWindow, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(mWindow, true);

        // Background Fill Color
        glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        update();
        render();

        // Flip Buffers and Draw
        glfwSwapBuffers(mWindow);
        glfwPollEvents();
    }  
    
    onExit();
    glfwTerminate();

    return EXIT_SUCCESS;
}

static void init()
{
    initSPH();
    initShaders();
    initCamera();
    initOfflineRendering();
    initPostprocessing();

    GLenum  err =  glGetError();
    if(err != GL_NO_ERROR)
    {
      std::cout << "OpenGL Error:" << err << std::endl;
    }

    return;
}

//
// Shader initializations
//

static void initShaders()
{
    init_map();
    init_pointSplat();
    init_gaussianBlur();
    init_depthview();
}

static void init_pointSplat()
{
  GLuint shaders[2];
  shaders[0] = ShaderUtils::loadShader(pointSplat_files[0], GL_VERTEX_SHADER);
  shaders[1] = ShaderUtils::loadShader(pointSplat_files[1], GL_FRAGMENT_SHADER);
  pointSplat_program = ShaderUtils::linkShaderProgram(&shaders[0], 2, true);

  //create vao
  glGenVertexArrays(1,&pointSplat_vao);
  glBindVertexArray(pointSplat_vao);

  //create vbo
  float* pos = NULL;
  int numVertex;
  simulator->getData(&pos, numVertex);   
  glGenBuffers(1,&vbo_pointSplat_pos);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_pointSplat_pos);
  glBufferStorage(GL_ARRAY_BUFFER, numVertex * 3 * sizeof(float), pos, GL_MAP_WRITE_BIT);
  
  //bind attributes
  glVertexAttribBinding(POINTSPLAT_LOCATION_POS, 0);
  glBindVertexBuffer(0, vbo_pointSplat_pos, 0, sizeof(float) * 3);
  glVertexAttribFormat(POINTSPLAT_LOCATION_POS, 3, GL_FLOAT, GL_FALSE, 0);
  glEnableVertexAttribArray(POINTSPLAT_LOCATION_POS);

  //recover
  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
 
} 

static void init_map()
{
  GLuint shaders[2];
  shaders[0] = ShaderUtils::loadShader(lighting_files[0], GL_VERTEX_SHADER);
  shaders[1] = ShaderUtils::loadShader(lighting_files[1], GL_FRAGMENT_SHADER);
  lighting_program = ShaderUtils::linkShaderProgram(&shaders[0], 2, true);

  //create vao
  glGenVertexArrays(1,&lighting_vao);
  glBindVertexArray(lighting_vao);

  //create vbo
  float* pos = NULL;
  int numVertex;
  simulator->getData(&pos, numVertex);   
  glGenBuffers(1,&vbo_lighting_pos);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_lighting_pos);
  glBufferStorage(GL_ARRAY_BUFFER, numVertex * 3 * sizeof(float), pos, GL_MAP_WRITE_BIT);
  
  //bind attributes
  glVertexAttribBinding(LIGHTING_LOCATION_POS, 0);
  glBindVertexBuffer(0, vbo_lighting_pos, 0, sizeof(float) * 3);
  glVertexAttribFormat(LIGHTING_LOCATION_POS, 3, GL_FLOAT, GL_FALSE, 0);
  glEnableVertexAttribArray(LIGHTING_LOCATION_POS);

  //recover
  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static void init_gaussianBlur()
{
  GLuint shaders[1];
  shaders[0] = ShaderUtils::loadShader(gaussianBlur_files[0], GL_COMPUTE_SHADER);  
  gaussianBlur_program = ShaderUtils::linkShaderProgram(shaders, 1, true);
}

static void init_depthview()
{
  GLuint shaders[2];
  shaders[0] = ShaderUtils::loadShader(depthview_files[0], GL_VERTEX_SHADER);
  shaders[1] = ShaderUtils::loadShader(depthview_files[1], GL_FRAGMENT_SHADER);
  depthview_program = ShaderUtils::linkShaderProgram(&shaders[0], 2, true);

  glGenVertexArrays(1,&depthview_vao);
  glBindVertexArray(depthview_vao);
  glGenBuffers(1,&vbo_depthview_pos);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_depthview_pos);


  GLfloat fullscreen[] = {-1.f, -1.f, 0,
                          -1.f, 1.f, 0,
                          1.f, -1.f, 0,
                          -1.f, 1.f, 0,
                          1.f, -1.f, 0,
                          1.f, 1.f, 0 };


  glBufferStorage(GL_ARRAY_BUFFER, sizeof(fullscreen), fullscreen, GL_MAP_READ_BIT);
  
  //bind attributes
  glVertexAttribBinding(DEPTHVIEW_LOCATION_POS, 0);
  glBindVertexBuffer(0, vbo_depthview_pos, 0, sizeof(float) * 3);
  glVertexAttribFormat(DEPTHVIEW_LOCATION_POS, 3, GL_FLOAT, GL_FALSE, 0);
  glEnableVertexAttribArray(DEPTHVIEW_LOCATION_POS);

  //recover
  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

}

//
//SPH Initialization
//

static void initSPH()
{
    delta = 1./60;

    box.spanX = 25;
    box.spanY = 25;
    box.spanZ = 25;

    SPHSim::SPHConfig conf;
    conf.kernel_h = KERNEL_H;
    conf.k = K;
    conf.density0 = RESET_DENSITY;
    conf.g = G;
    conf.miu = MIU;
    conf.boundarydamping = BOUNDARYDAMPING;
    conf.damping = DAMPING;
    conf.tensionCoe = TENSION_COEF;
    conf.box = box;
    simulator = new SPHSim::SPHSimulator(conf);
    simulator->setup();

    return;
}

//
//Camera Initialization
//

static void initCamera()
{
    m2w = rotate(0.5f, vec3(1,0,0));
    m2w = rotate(0.1f, vec3(0,0,1));

    //set up camera
    vec3 eye(0,30,70);
    vec3 center(0,0,0);
    vec3 up(0,1,0);
    w2v = lookAt(eye, center, up);

    //set up perspective projection
    pers_proj = glm::perspective(45.0f, ((float)mWidth)/mHeight, 0.1f, 500.0f);    

    //setup cursor position
    cursor_cur_x = CURSOR_POS_INVALID;
    cursor_cur_y = CURSOR_POS_INVALID;
    cursor_pre_x = CURSOR_POS_INVALID;
    cursor_pre_y = CURSOR_POS_INVALID;
    
    return;
}

//
// Offline Rendering Initialization: fbo
//

static void initOfflineRendering()
{

    //create fbo
    glCreateFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    //create color buffer texture
    glGenTextures(1, &normalBuffer);
    glBindTexture(GL_TEXTURE_2D, normalBuffer);
    glTextureStorage2D(normalBuffer, 1, GL_RGBA32F, mWidth, mHeight);
    
    //attach the color buffer texture to fbo
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, normalBuffer, 0);

    //create depth buffer texture
    glGenTextures(1, &depthBuffer);
    glBindTexture(GL_TEXTURE_2D, depthBuffer);
    glTextureStorage2D(depthBuffer, 1, GL_DEPTH_COMPONENT32F, mWidth, mHeight);

    //attach buffer textures to fbo
    glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT , depthBuffer, 0);
    
    //tell opengl which color attachments to render to
    GLuint buffers[] = {GL_COLOR_ATTACHMENT0};
    glDrawBuffers(1, buffers);

    //recover
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    //check fbo completeness
    GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
    if(status!=GL_FRAMEBUFFER_COMPLETE)
    {
       printf("bad!!\n");
    }

    return;
}


//
// initialization for post processing
// eg. textures, buffers, ect.
//
static void initPostprocessing()
{
    glGenTextures(1, &blurredNormal);
    glBindTexture(GL_TEXTURE_2D, blurredNormal);
    glTextureStorage2D(blurredNormal, 1, GL_RGBA32F, mWidth, mHeight);  
}

//
// RENDERING
//
static void render()
{
    renderBox();

    Timer timer;
    timer.start();
    //renderParticlesAsPoints();
    renderParticles(true);
    renderFluid();
    timer.stop();
    double time = timer.duration();
    std::cout << "render: " << time * 1000 << "ms" << std::endl;

}

static void renderBox()
{
    setTransform();

    double boxX = box.spanX;
    double boxY = box.spanY;
    double boxZ = box.spanZ;

  glBegin(GL_LINES);
      glColor3f(1,0,0);
      glVertex3f( -boxX/2, -boxY/2,  boxZ/2);
      glVertex3f(  boxX/2, -boxY/2,  boxZ/2);
      glVertex3f( -boxX/2,  boxY/2,  boxZ/2);
      glVertex3f(  boxX/2,  boxY/2,  boxZ/2);
      glVertex3f( -boxX/2, -boxY/2,  -boxZ/2);
      glVertex3f(  boxX/2, -boxY/2,  -boxZ/2);
      glVertex3f( -boxX/2,  boxY/2,  -boxZ/2);
      glVertex3f(  boxX/2,  boxY/2,  -boxZ/2);

      glVertex3f(  boxX/2,  boxY/2,   boxZ/2);
      glVertex3f(  boxX/2, -boxY/2,   boxZ/2);
      glVertex3f(  boxX/2,  boxY/2,  -boxZ/2);
      glVertex3f(  boxX/2, -boxY/2,  -boxZ/2);
      glVertex3f( -boxX/2,  boxY/2,   boxZ/2);
      glVertex3f( -boxX/2, -boxY/2,   boxZ/2);
      glVertex3f( -boxX/2,  boxY/2,  -boxZ/2);
      glVertex3f( -boxX/2, -boxY/2,  -boxZ/2);

      glVertex3f(  boxX/2,  boxY/2,   boxZ/2);
      glVertex3f(  boxX/2,  boxY/2,  -boxZ/2);
      glVertex3f(  boxX/2,  -boxY/2,  boxZ/2);
      glVertex3f(  boxX/2,  -boxY/2,  -boxZ/2);
      glVertex3f(  -boxX/2,  boxY/2,  boxZ/2);
      glVertex3f(  -boxX/2,  boxY/2,  -boxZ/2);
      glVertex3f(  -boxX/2,  -boxY/2, boxZ/2);
      glVertex3f(  -boxX/2,  -boxY/2,  -boxZ/2);

  glEnd();
}

static void renderParticles(bool offscreen)
{
    if(offscreen)
    {
      //use fbo
      glBindFramebuffer(GL_FRAMEBUFFER, fbo);
      glViewport(0,0,mWidth, mHeight);
      glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glEnable(GL_DEPTH_TEST);
    }
    //use program
    glUseProgram(pointSplat_program);
    glBindVertexArray(pointSplat_vao);
    
    //set uniforms
    glUniform1f(POINTSPLAT_LOCATION_R, PARTICLE_SPLAT_SIZE);
    glUniformMatrix4fv(POINTSPLAT_LOCATION_M2W,        1, GL_FALSE, value_ptr(m2w));
    glUniformMatrix4fv(POINTSPLAT_LOCATION_W2V, 	      1, GL_FALSE, value_ptr(w2v));
    glUniformMatrix4fv(POINTSPLAT_LOCATION_PERSP_PROJ, 1, GL_FALSE, value_ptr(pers_proj));

    //set configuration
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_POINT_SPRITE);

    //update position vbo
    int numVertex;
    float* buffer;
    buffer = (float*) glMapNamedBuffer(vbo_pointSplat_pos,  GL_WRITE_ONLY);
    simulator->getData(&buffer, numVertex);   
    glUnmapNamedBuffer(vbo_pointSplat_pos);

    //draw
    glDrawArrays(GL_POINTS, 0, numVertex);

    //reset
    glDisable(GL_PROGRAM_POINT_SIZE);
    glDisable(GL_POINT_SPRITE);
    glUseProgram(0);
    glBindVertexArray(0);

    if(offscreen)
    {
      glViewport(0,0, mWidth, mHeight);
      glBindFramebuffer(GL_FRAMEBUFFER, 0);  

      for(int i=0; i<4; i++)
      {
        blurTexture(normalBuffer, blurredNormal, mWidth, mHeight);
        blurTexture(blurredNormal, normalBuffer, mWidth, mHeight);
      }
      blurTexture(normalBuffer, blurredNormal, mWidth, mHeight);
    }
}

static void renderParticlesAsPoints()
{
     float* data = NULL;
     int num;
     simulator->getData(&data, num);
 
     glPointSize(1);
     glBegin(GL_POINTS);
     glColor3f(1,0,0);
     for(int i=0; i<num; i++)
     {
         glVertex3f(data[3*i + 0],
                    data[3*i + 1],
                    data[3*i + 2]);
     }
     glEnd();
     free(data);
}

static void renderFluid()
{
    //use program
    glUseProgram(lighting_program);
    glBindVertexArray(lighting_vao);
    
    //set uniforms
    glUniform1f(LIGHTING_LOCATION_R, PARTICLE_SPLAT_SIZE);
    glUniformMatrix4fv(LIGHTING_LOCATION_M2W,        1, GL_FALSE, value_ptr(m2w));
    glUniformMatrix4fv(LIGHTING_LOCATION_W2V, 	    1, GL_FALSE, value_ptr(w2v));
    glUniformMatrix4fv(LIGHTING_LOCATION_PERSP_PROJ, 1, GL_FALSE, value_ptr(pers_proj));
    glUniform3fv(LIGHTING_LOCATION_AMBIENT, 1, value_ptr(ambient));
    glUniform3fv(LIGHTING_LOCATION_DIFFUSE, 1, value_ptr(diffuse));
    glUniform3fv(LIGHTING_LOCATION_SPECULAR, 1, value_ptr(specular));
    vec3 lDirV = vec3(w2v * vec4(lightDir,0));
    glUniform3fv(LIGHTING_LOCATION_LDIR, 1, value_ptr(lDirV));
    glUniform1f(LIGHTING_LOCATION_SHINESNESS, shineness);

    //set configuration
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_POINT_SPRITE);

    //update position vbo
    int numVertex;
    float* buffer;
    buffer = (float*) glMapNamedBuffer(vbo_lighting_pos,  GL_WRITE_ONLY);
    simulator->getData(&buffer, numVertex);
    glUnmapNamedBuffer(vbo_lighting_pos);
   
    //bind texture
    glBindTexture(GL_TEXTURE_2D, blurredNormal);

    //draw
    glDrawArrays(GL_POINTS, 0, numVertex);

    //reset
    glDisable(GL_PROGRAM_POINT_SIZE);
    glDisable(GL_POINT_SPRITE);
    glUseProgram(0);
    glBindVertexArray(0);
}

static void blurTexture(GLuint input, GLuint output, int width, int height)
{
   glUseProgram(gaussianBlur_program);
   glUniform1i(0, width);
   glUniform1i(1, height);
   glBindImageTexture(0, input , 0, GL_FALSE, 0, GL_READ_ONLY, GL_RGBA32F);
   glBindImageTexture(1, output, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
   glDispatchCompute((width+31)/32,(height+31)/32,1);
   glUseProgram(0);

   return;
}

static void setTransform()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixf(value_ptr(w2v));
    glMultMatrixf(value_ptr(m2w));

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMultMatrixf(value_ptr(pers_proj));

    glMatrixMode(GL_MODELVIEW);
}

//
// GAME LOGIC UPDATES
//
static void update()
{
    if(is_mouse_pressed)
        duringPress();

    vec4 g(0.f,-1.f,0.f,0.f);
    g = inverse(m2w) * g;
    simulator->setGravityDirection(g.x, g.y, g.z);

    Timer timer;
    timer.start();
    simulator->update(delta);
    timer.stop(); 
    double time = timer.duration();
    std::cout << "simulation:" << time*1000 << "ms" << std::endl;
}

//
// MOUSE INPUT HANDLING
//

static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        if(action == GLFW_PRESS)
            onPress();
        else if(action == GLFW_RELEASE)
            onRelease();
    }
}

static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
   w2v = translate(w2v,vec3(0,0,yoffset));
}

static void onPress()
{
    is_mouse_pressed = true;
}

static void onRelease()
{
    is_mouse_pressed = false;
    cursor_cur_x = CURSOR_POS_INVALID;
    cursor_cur_y = CURSOR_POS_INVALID;
    cursor_pre_x = CURSOR_POS_INVALID;
    cursor_pre_y = CURSOR_POS_INVALID;
}

static void duringPress()
{
    //update cursor position
    cursor_pre_x = cursor_cur_x; cursor_pre_y = cursor_cur_y;
    glfwGetCursorPos(mWindow, &cursor_cur_x, &cursor_cur_y);
    if(cursor_pre_x == CURSOR_POS_INVALID)
    {
        cursor_pre_x = cursor_cur_x;
        cursor_pre_y = cursor_cur_y;
    }

    //arcball update
    arcball_update();
}


//
// ARCBALL MECHANISM
//
// ref: http://courses.cms.caltech.edu/cs171/assignments/hw3/hw3-notes/notes-hw3.html
//
static void arcball_update()
{
    //not need to update if mouse is not moved
    if(cursor_pre_x == cursor_cur_x && cursor_pre_y == cursor_cur_y)
        return;

    //compute the position in NDC space
    double ndc_pre_x = (cursor_pre_x - mWidth/2) / (mWidth/2);
    double ndc_pre_y = (mHeight/2 -cursor_pre_y) / (mHeight/2);
    double ndc_cur_x = (cursor_cur_x - mWidth/2) / (mWidth/2);
    double ndc_cur_y = (mHeight/2 -cursor_cur_y) / (mHeight/2);

    //estimate the z coordinates in NDC space using the sphere surface constraints
    double tmp, ndc_pre_z, ndc_cur_z;

    tmp = ndc_pre_x * ndc_pre_x + ndc_pre_y * ndc_pre_y;
    if(tmp > 1)
        ndc_pre_z = 0;
    else
        ndc_pre_z = std::sqrt(1 - tmp);

    tmp = ndc_cur_x * ndc_cur_x + ndc_cur_y * ndc_cur_y;
    if(tmp > 1)
        ndc_cur_z = 0;
    else
        ndc_cur_z = std::sqrt(1 - tmp);

    //compute vector u and rotate angle
    vec3 pre(ndc_pre_x, ndc_pre_y, ndc_pre_z);
    vec3 cur(ndc_cur_x, ndc_cur_y, ndc_cur_z);
    vec3 u = normalize(cross(pre, cur));
    float r = std::acos(std::min(dot(pre,cur),1.f));

    //multiplies speed factor
    r *= speed_factor;

    //update m2w matrix
    m2w = rotate(r, u) * m2w;

    return;
}

static void onExit()
{
    glDeleteProgram(pointSplat_program);
    glDeleteVertexArrays(1,&pointSplat_vao);
    glDeleteBuffers(1, &vbo_pointSplat_pos);
    delete simulator;
}

static void outputTexture2File(GLuint texture, const char* file, GLenum format, GLenum type, int channel)
{ 
   unsigned char* data = new unsigned char[mWidth * mHeight * channel];

   glBindTexture(GL_TEXTURE_2D, texture);

   //get texture image
   glGetTexImage(GL_TEXTURE_2D, 0,
  	format,
  	type,
  	data);
   
   glBindTexture(GL_TEXTURE_2D, 0);
   
   //write data to image file
   int save_result = SOIL_save_image
    (
        file,
        SOIL_SAVE_TYPE_PNG,
        mWidth, mHeight, channel,
        data
    );

    delete[] data;
   
   return;
}

static void visualizeDepth(GLuint depth)
{
  glUseProgram(depthview_program);
  glBindVertexArray(depthview_vao);
  //glBindImageTexture(0, depthBuffer, 0, GL_FALSE, 0, GL_READ_ONLY, GL_R32F);
  glBindTexture(GL_TEXTURE_2D, depthBuffer);
  glDrawArrays(GL_TRIANGLES, 0, 6);
  glUseProgram(0);
}


static void testLoadImage()
{
     /* load an image file directly as a new OpenGL texture */
     GLuint tex_2d = SOIL_load_OGL_texture
     (
      "img.png",
      SOIL_LOAD_AUTO,
      SOIL_CREATE_NEW_ID,
      SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT
     );

     /* check for an error during the load process */
    if( 0 == tex_2d )
    {
      printf( "SOIL loading error: '%s'\n", SOIL_last_result() );
    } 
}
