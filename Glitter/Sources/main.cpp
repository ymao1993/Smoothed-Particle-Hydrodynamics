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

//Rendering
static const char* shaderFiles[] =
{
 "../Glitter/Shaders/pointSplat.vs",
 "../Glitter/Shaders/pointSplat.fs",
 "../Glitter/Shaders/GaussianBlur.cs"
};
static GLuint program;
static GLuint vao;
static GLuint vbo_pos;
static GLuint colorBuffer;
static GLuint depthBuffer;
static GLuint fbo;
static GLuint blurredColor;

//shaders
static GLuint blur_program;

//forward declarations
static void init();
static void render();
static void update();
static void onExit();

static void initShaders();
static void initSPH();
static void initCamera();
static void initOfflineRendering();
static void initPostprocessing();
static void renderBox();
static void renderParticles(bool offscreen);
static void blurTexture(GLuint input, GLuint output, int width, int height);
static void setTransform();
static void outputTexture2File(GLuint texture, const char* file);

//tests
static void testLoadImage();

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

    return;
}

enum
{
  LOCATION_POS,
};

static void initShaders()
{   
    //build shader program
    GLuint shaders[2];
    shaders[0] = ShaderUtils::loadShader(shaderFiles[0], GL_VERTEX_SHADER);
    shaders[1] = ShaderUtils::loadShader(shaderFiles[1], GL_FRAGMENT_SHADER);
    shaders[2] = ShaderUtils::loadShader(shaderFiles[2], GL_COMPUTE_SHADER);
    program = ShaderUtils::linkShaderProgram(shaders, 2, true);

    blur_program = ShaderUtils::linkShaderProgram(&shaders[2], 1, true);

    //create vao
    glGenVertexArrays(1,&vao);
    glBindVertexArray(vao);

    //create vbo
    float* pos = NULL;
    int numVertex;
    simulator->getData(&pos, numVertex);   
    glGenBuffers(1,&vbo_pos);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_pos);
    glBufferStorage(GL_ARRAY_BUFFER, numVertex * 3 * sizeof(float), pos, GL_MAP_WRITE_BIT);
    
    //bind attributes
    glVertexAttribBinding(LOCATION_POS, 0);
    glBindVertexBuffer(0, vbo_pos, 0, sizeof(float) * 3);
    glVertexAttribFormat(LOCATION_POS, 3, GL_FLOAT, GL_FALSE, 0);
    glEnableVertexAttribArray(LOCATION_POS);

}
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
static void initCamera()
{
    //set up camera
    vec3 eye(0,0,30);
    vec3 center(0,0,0);
    vec3 up(0,1,0);
    w2v = lookAt(eye, center, up);

    //set up perspective projection
    pers_proj = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);    

    //setup cursor position
    cursor_cur_x = CURSOR_POS_INVALID;
    cursor_cur_y = CURSOR_POS_INVALID;
    cursor_pre_x = CURSOR_POS_INVALID;
    cursor_pre_y = CURSOR_POS_INVALID;
    
    return;
}

static void initOfflineRendering()
{

    //create fbo
    glCreateFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    //create color buffer texture
    glGenTextures(1, &colorBuffer);
    glBindTexture(GL_TEXTURE_2D, colorBuffer);
    glTextureStorage2D(colorBuffer, 1, GL_RGBA8, 512, 512);
    
    //attach the color buffer texture to fbo
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, colorBuffer, 0);

    //create depth buffer texture
    glGenTextures(1, &depthBuffer);
    glBindTexture(GL_TEXTURE_2D, depthBuffer);
    glTextureStorage2D(depthBuffer, 1, GL_DEPTH_COMPONENT32F, 512, 512);

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

static void initPostprocessing()
{
    glGenTextures(1, &blurredColor);
    glBindTexture(GL_TEXTURE_2D, blurredColor);
    glTextureStorage2D(blurredColor, 1, GL_RGBA8, 512, 512);   
}

//
// RENDERING
//
static void render()
{
    renderBox();
    renderParticles(false);
    renderParticles(true);
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

enum 
{
   LOCATION_R,
   LOCATION_M2W,
   LOCATION_W2V,
   LOCATION_PERSP_PROJ,
};


static void renderParticles(bool offscreen)
{
    if(offscreen)
    {
      //use fbo
      glBindFramebuffer(GL_FRAMEBUFFER, fbo);
      glViewport(0,0,512,512);
      glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
    //use program
    glUseProgram(program);
    glBindVertexArray(vao);
    
    //set uniforms
    glUniform1f(LOCATION_R, 500);
    glUniformMatrix4fv(LOCATION_M2W,        1, GL_FALSE, value_ptr(m2w));
    glUniformMatrix4fv(LOCATION_W2V, 	    1, GL_FALSE, value_ptr(w2v));
    glUniformMatrix4fv(LOCATION_PERSP_PROJ, 1, GL_FALSE, value_ptr(pers_proj));

    //set configuration
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_POINT_SPRITE);

    //update position vbo
    int numVertex;
    float* buffer;
    buffer = (float*) glMapBuffer(GL_ARRAY_BUFFER,  GL_WRITE_ONLY);
    simulator->getData(&buffer, numVertex);   
    glUnmapBuffer(GL_ARRAY_BUFFER);

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

      blurTexture(colorBuffer, blurredColor, 512, 512);

      outputTexture2File(colorBuffer, "color.png");
      outputTexture2File(blurredColor, "bcolor.png");  
    }
}

static void blurTexture(GLuint input, GLuint output, int width, int height)
{
  
   glUseProgram(blur_program);   
   glUniform1i(0, width);
   glUniform1i(1, height);

   glBindImageTexture(0, input , 0, GL_FALSE, 0, GL_READ_ONLY, GL_RGBA8);
   glBindImageTexture(1, output, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);
   glDispatchCompute(width/32,height/32,1);
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
    simulator->update(delta);
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
    glDeleteProgram(program);
    glDeleteVertexArrays(1,&vao);
    glDeleteBuffers(1, &vbo_pos);
    delete simulator;
}

static void outputTexture2File(GLuint texture, const char* file)
{ 
   static unsigned char* data = new unsigned char[512 * 512 * 4];

   glBindTexture(GL_TEXTURE_2D, texture);

   //get texture image
   glGetTexImage(GL_TEXTURE_2D, 0,
  	GL_RGBA,
  	GL_UNSIGNED_BYTE,
  	data);
   
   glBindTexture(GL_TEXTURE_2D, 0);
   
   //write data to image file
   int save_result = SOIL_save_image
    (
        file,
        SOIL_SAVE_TYPE_PNG,
        512, 512, 4,
        data
    );
   
   return;
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
