#include "FIELD_2D.h"
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#elif __linux__
#include <GL/glut.h>
#endif

using namespace std;

// --- EXTREME BRANCHING SETTINGS ---
float length = 40.0f;      
float eps_bar = 0.4f;     // thinner interface, less stabilization
float dee = 0.4f;        // stronger anisotropy — tips WANT to form
int j = 6;
float tau = 0.3f;         // scale with eps_bar^2 roughly
float lambda_param = 4.0f; // moderate driving
float L_sat = 2.0f;
float dt = 0.003f;

#pragma region OpenGL Init
// Text for the title bar of the window
string windowLabel("Phase Field Model");

// mouse tracking variables
int xMouse = -1;
int yMouse = -1;
int mouseButton = -1;
int mouseState = -1;
int mouseModifiers = -1;

// current grid cell the mouse is pointing at
int xField = -1;
int yField = -1;

// animate the current runEverytime()?
bool animate = true;

// the resolution of the OpenGL window -- independent of the field resolution
int xScreenRes = 750;
int yScreenRes = 750;

// draw the grid over the field?
bool drawingGrid = false;

// print out what the mouse is pointing at?
bool drawingValues = true;

// currently capturing frames for a movie?
bool captureMovie = false;

// the current viewer eye position
float eyeCenter[] = {0.5, 0.5, 1};

// current zoom level into the field
float zoom = 1.0;


#pragma endregion

// resolution of the field
int xRes = 200,    // horizontal in glass
    zRes = 200,    // vertical in glass
    yLayers = 30;  // number of slabs in 3d

std::vector<FIELD_2D> vapor_u; // stack of 2D layers to represent 3D field
FIELD_2D phaseField(xRes, zRes); // the 2d substrate phase field, livse only at layer 0.
FIELD_2D KBYS_dpdy(xRes, zRes);
FIELD_2D KBYS_dpdx(xRes, zRes);
FIELD_2D heightField(xRes, zRes);
float heightScale = 1.0f; // visualization scaling, tweak

float delta = length / (xRes - 1);

// the field being drawn and manipulated
FIELD_2D field(xRes, zRes);

float *aData, *bData, *newAData, *newBData;

void runOnce();
void runEverytime();

///////////////////////////////////////////////////////////////////////
// Figure out which field element is being pointed at, set xField and
// yField to them 
///////////////////////////////////////////////////////////////////////
void refreshMouseFieldIndex(int x, int y)
{
  // make the lower left the origin
  y = yScreenRes - y;

  float xNorm = (float)x / xScreenRes;
  float yNorm = (float)y / yScreenRes;

  float halfZoom = 0.5 * zoom;
  float xWorldMin = eyeCenter[0] - halfZoom;
  float xWorldMax = eyeCenter[0] + halfZoom;

  // get the bounds of the field in screen coordinates
  //
  // if non-square textures are ever supported, change the 0.0 and 1.0 below
  float xMin = (0.0 - xWorldMin) / (xWorldMax - xWorldMin);
  float xMax = (1.0 - xWorldMin) / (xWorldMax - xWorldMin);

  float yWorldMin = eyeCenter[1] - halfZoom;
  float yWorldMax = eyeCenter[1] + halfZoom;

  float yMin = (0.0 - yWorldMin) / (yWorldMax - yWorldMin);
  float yMax = (1.0 - yWorldMin) / (yWorldMax - yWorldMin);

  float xScale = 1.0;
  float yScale = 1.0;

  if (xRes < zRes)
    xScale = (float)zRes / xRes;
  if (xRes > zRes)
    yScale = (float)xRes / zRes;

  // index into the field after normalizing according to screen
  // coordinates
  xField = xScale * xRes * ((xNorm - xMin) / (xMax - xMin));
  yField = yScale * zRes * ((yNorm - yMin) / (yMax - yMin));

  // clamp to something inside the field
  xField = (xField < 0) ? 0 : xField;
  xField = (xField >= xRes) ? xRes - 1 : xField;
  yField = (yField < 0) ? 0 : yField;
  yField = (yField >= zRes) ? zRes - 1 : yField;
}

///////////////////////////////////////////////////////////////////////
// Print a string to the GL window
///////////////////////////////////////////////////////////////////////
void printGlString(string output)
{
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
  for (unsigned int x = 0; x < output.size(); x++)
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, output[x]);
}

///////////////////////////////////////////////////////////////////////
// dump the field contents to a GL texture for drawing
///////////////////////////////////////////////////////////////////////
void updateTexture(FIELD_2D &texture)
{
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, texture.xRes(), texture.zRes(), 0, GL_LUMINANCE, GL_FLOAT, texture.data());

  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glEnable(GL_TEXTURE_2D);
}

///////////////////////////////////////////////////////////////////////
// draw a grid over everything
///////////////////////////////////////////////////////////////////////
void drawGrid()
{
  glColor4f(0.1, 0.1, 0.1, 1.0);

  float dx = 1.0 / xRes;
  float dy = 1.0 / zRes;

  if (xRes < zRes)
    dx *= (float)xRes / zRes;
  if (xRes > zRes)
    dy *= (float)zRes / xRes;

  glBegin(GL_LINES);
  for (int x = 0; x < field.xRes() + 1; x++) {
    glVertex3f(x * dx, 0, 1);
    glVertex3f(x * dx, 1, 1);
  }
  for (int y = 0; y < field.zRes() + 1; y++) {
    glVertex3f(0, y * dy, 1);
    glVertex3f(1, y * dy, 1);
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  // Make ensuing transforms affect the projection matrix
  glMatrixMode(GL_PROJECTION);

  // set the projection matrix to an orthographic view
  glLoadIdentity();
  float halfZoom = zoom * 0.5;

  glOrtho(-halfZoom, halfZoom, -halfZoom, halfZoom, -10, 10);

  // set the matrix mode back to modelview
  glMatrixMode(GL_MODELVIEW);

  // set the lookat transform
  glLoadIdentity();
  gluLookAt(eyeCenter[0], eyeCenter[1], 1, // eye
            eyeCenter[0], eyeCenter[1], 0, // center
            0, 1, 0);                      // up

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  float xLength = 1.0;
  float yLength = 1.0;

  if (xRes < zRes)
    xLength = (float)xRes / zRes;
  if (zRes < xRes)
    yLength = (float)zRes / xRes;

  glEnable(GL_TEXTURE_2D);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glTexCoord2f(0.0, 1.0);
  glVertex3f(0.0, yLength, 0.0);
  glTexCoord2f(1.0, 1.0);
  glVertex3f(xLength, yLength, 0.0);
  glTexCoord2f(1.0, 0.0);
  glVertex3f(xLength, 0.0, 0.0);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  // draw the grid, but only if the user wants
  if (drawingGrid)
    drawGrid();

  // if there's a valid field index, print it
  if (xField >= 0 && yField >= 0 && xField < field.xRes() && yField < field.zRes()) {
    glLoadIdentity();

    // must set color before setting raster position, otherwise it won't take
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

    // normalized screen coordinates (-0.5, 0.5), due to the glLoadIdentity
    float halfZoom = 0.5 * zoom;
    glRasterPos3f(-halfZoom * 0.95, -halfZoom * 0.95, 0);

    // build the field value string
    char buffer[256];
    string fieldValue("(");
    sprintf(buffer, "%i", xField);
    fieldValue = fieldValue + string(buffer);
    sprintf(buffer, "%i", yField);
    fieldValue = fieldValue + string(", ") + string(buffer) + string(") = ");
    sprintf(buffer, "%f", field(xField, yField));
    fieldValue = fieldValue + string(buffer);

    // draw the grid, but only if the user wants
    if (drawingValues)
      printGlString(fieldValue);
  }

  glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////
void printCommands()
{
  cout << "=============================================================== " << endl;
  cout << " Field viewer code CPSC 679" << endl;
  cout << "=============================================================== " << endl;
  cout << " q           - quit" << endl;
  cout << " v           - type the value of the cell under the mouse" << endl;
  cout << " g           - throw a grid over everything" << endl;
  cout << " w           - write out a PPM file " << endl;
  cout << " left mouse  - pan around" << endl;
  cout << " right mouse - zoom in and out " << endl;
  cout << " shift left mouse - draw on the grid " << endl;
}

///////////////////////////////////////////////////////////////////////
// Map the arrow keys to something here
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y)
{
  switch (key) {
  case GLUT_KEY_LEFT:
    break;
  case GLUT_KEY_RIGHT:
    break;
  case GLUT_KEY_UP:
    break;
  case GLUT_KEY_DOWN:
    break;
  default:
    break;
  }
}

///////////////////////////////////////////////////////////////////////
// Map the keyboard keys to something here
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'a':
    animate = !animate;
    break;
  case 'g':
    drawingGrid = !drawingGrid;
    break;
  case '?':
    printCommands();
    break;
  case 'v':
    drawingValues = !drawingValues;
    break;
  case 'w': {
    static int count = 0;
    char buffer[256];
    sprintf(buffer, "output_%i.ppm", count);
    field.writePPM(buffer);
    count++;
  } break;
  case 'q':
    exit(0);
    break;
  default:
    break;
  }
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is clicked
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y)
{
  int modifiers = glutGetModifiers();
  mouseButton = button;
  mouseState = state;
  mouseModifiers = modifiers;

  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && modifiers & GLUT_ACTIVE_SHIFT) {
    // figure out which cell we're pointing at
    refreshMouseFieldIndex(x, y);

    // set the cell
    for (int y = yField - 4; y < yField + 5; y++) {
      for (int x = xField - 4; x < xField + 5; x++) {
        // bField(x, y) = 0;
      }
    }

    // make sure nothing else is called
    return;
  }

  xMouse = x;
  yMouse = y;
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is clicked and moving
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y)
{
  if (mouseButton == GLUT_LEFT_BUTTON && mouseState == GLUT_DOWN && mouseModifiers & GLUT_ACTIVE_SHIFT) {
    // figure out which cell we're pointing at
    refreshMouseFieldIndex(x, y);

    // set the cell
    // for (int y = yField - 4; y < yField + 5; y++) {
    //   for (int x = xField - 4; x < xField + 5; x++) {
    //     bField(x, y) = 0;
    //   }
    // }

    // make sure nothing else is called
    return;
  }

  float xDiff = x - xMouse;
  float yDiff = y - yMouse;
  float speed = 0.001;

  if (mouseButton == GLUT_LEFT_BUTTON) {
    eyeCenter[0] -= xDiff * speed;
    eyeCenter[1] += yDiff * speed;
  }
  if (mouseButton == GLUT_RIGHT_BUTTON)
    zoom -= yDiff * speed;

  xMouse = x;
  yMouse = y;
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is not clicked and moving
///////////////////////////////////////////////////////////////////////
void glutPassiveMouseMotion(int x, int y)
{
  refreshMouseFieldIndex(x, y);
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate) {
    runEverytime();
  }
  updateTexture(field);
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowSize(xScreenRes, yScreenRes);
  glutInitWindowPosition(10, 10);
  glutCreateWindow(windowLabel.c_str());

  // set the viewport resolution (w x h)
  glViewport(0, 0, (GLsizei)xScreenRes, (GLsizei)yScreenRes);

  // set the background color to gray
  glClearColor(0.1, 0.1, 0.1, 0);

  // register all the callbacks
  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);
  glutPassiveMotionFunc(&glutPassiveMouseMotion);

  // enter the infinite GL loop
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  // In case the field is rectangular, make sure to center the eye
  if (xRes < zRes) {
    float xLength = (float)xRes / zRes;
    eyeCenter[0] = xLength * 0.5;
  }
  if (zRes < xRes) {
    float yLength = (float)zRes / xRes;
    eyeCenter[1] = yLength * 0.5;
  }
  
  printCommands();

  cout << "Simulating " << endl;

  runOnce();

  // initialize GLUT and GL
  glutInit(&argc, argv);

  // open the GL window
  glvuWindow();
  return 1;
}

///////////////////////////////////////////////////////////////////////
// This is called once at the beginning so you can precache
// something here
///////////////////////////////////////////////////////////////////////
void runOnce()
{
  for (int y = 0; y < yLayers; y++) {
    vapor_u.push_back(FIELD_2D(xRes, zRes));
  }

  heightField = 0.0f;
  phaseField = -1.0f; // everything is vapor initially

  int centerX = xRes / 2;
  int centerZ = zRes / 2;
  float radius = 1.0f;
  
  for (int i = 0; i < xRes; i++) {
    for (int k = 0; k < zRes; k++) {
      float dist = sqrt(pow(i - centerX, 2) + pow(k - centerZ, 2));
      // Smoothly ramps from +1 inside to -1 outside across the boundary
      phaseField(i, k) = -tanh((dist - radius) / (sqrt(2.0f) * eps_bar)); 
    }
  }

  // vapor gradient: glass is cold, air high above has lots of humidity
  for (int y = 0; y < yLayers; y++) {
    float height_fraction = (float)y / (float)(yLayers-1);

    for (int i = 0; i < xRes; i++) 
      for(int k = 0; k < zRes; k++) {
        vapor_u[y](i,k) = 1.0f;
      }
  }
}

float get_q(int x, int z, int y) {
  // Block diffusion through ice in the first few layers
  int reactive_layers = 3; 
  if (y < reactive_layers) {
    float phi = phaseField(x, z);
    return (1.0f - phi) / 2.0f; // 0 in solid, 1 in vapor
  }
  return 1.0f; // free air above
}

void create_toroid(FIELD_2D &f) {
  // left and right border
  for (int i = 1; i < zRes - 2; i++) {
    f(0,i) = f(xRes-2, i);
    f(xRes-1,i) = f(1, i);
  }
  // bottom and top border
  for (int i = 1; i < xRes - 2; i++) {
    f(i,0) = f(i, zRes-2);
    f(i,zRes-1) = f(i, 1);
  }
}

float space_derivative(const FIELD_2D &f, int dim, int x, int y, float delta) {
  // dim = 0 -> wrt x
  if (!dim) {
    return (f(x+1,y) - f(x-1,y))/(2*delta);
  }
  // dim = 1 -> wrt y
  else {
    return (f(x,y+1) - f(x,y-1))/(2*delta);
  }
}

// from Kim 2003
float find_theta(FIELD_2D &f, int x, int y, float delta) {
  float dp_dx = space_derivative(f, 0, x, y, delta);
  float dp_dy = space_derivative(f, 1, x, y, delta);

  float theta = -1.0 * (atan2(dp_dx, dp_dy) + M_PI/2);

  return theta;
}

void fill_KBYS_field(FIELD_2D &f, FIELD_2D &pField, int dim, float delta) {
  float dp, theta, KBYS_eps, de_dtheta, result;

  for (int y = 1; y < zRes-1; y++) {
    for (int x = 1; x < xRes-1; x++) {
      dp = space_derivative(pField, dim, x, y, delta);

      theta = find_theta(pField, x, y, delta); // in radians
      KBYS_eps = eps_bar * (1 + dee * cos(j * theta));
      de_dtheta = -eps_bar * dee * j * sin(j * theta);

      result = KBYS_eps * de_dtheta * dp;

      f(x,y) = result;
    }
  }

  create_toroid(f);
}

void runEverytime()
{
  std::vector<FIELD_2D> next_u = vapor_u; 
  FIELD_2D next_phase = phaseField;
  FIELD_2D dphi_dt_field(xRes, zRes); // To store the growth rate
  
  float dx2 = delta * delta;

  // calculate the stable anisotropic torque fields 
  fill_KBYS_field(KBYS_dpdy, phaseField, 1, delta); 
  fill_KBYS_field(KBYS_dpdx, phaseField, 0, delta);

  // STEP 1: evolve 2D Phase Field 
  for (int x = 1; x < xRes - 1; x++) {
    for (int z = 1; z < zRes - 1; z++) {
      
      float phi = phaseField(x, z);
      float u_s = vapor_u[0](x, z); 
      
      // get the stable Kobayashi anisotropic Laplacian
      float theta = find_theta(phaseField, x, z, delta);
      float KBYS_eps = eps_bar * (1.0f + dee * cos(j * theta)); // j=6 for ice
      
      float d_dx = space_derivative(KBYS_dpdy, 0, x, z, delta);
      float d_dz = space_derivative(KBYS_dpdx, 1, x, z, delta);   
      
      float laplace = (phaseField(x+1, z) + phaseField(x-1, z) + 
                       phaseField(x, z+1) + phaseField(x, z-1) - 
                       4.0f * phi) / dx2;

      // anisotropic curvature term
      float stable_anisotropy = (KBYS_eps * KBYS_eps * laplace) + d_dx - d_dz;

      // double well + tnterface coupling
      float f_prime = -phi + (phi * phi * phi);
      float g_prime = pow(1.0f - (phi * phi), 2);
      
      float random_val = ((rand() % 1000) / 1000.0f - 0.5f);
      float noise = random_val * 0.0f * g_prime; 
      
      float dphi_dt = (stable_anisotropy - f_prime + (lambda_param + noise) * g_prime * u_s) / tau;

      float grow = std::max(0.0f, dphi_dt);         // only accumulate solidification
      float mask = g_prime;                         // interface-only
      heightField(x,z) += dt * heightScale * grow * mask;
            
      dphi_dt_field(x, z) = dphi_dt;
      next_phase(x, z) = phi + dt * dphi_dt;
    }
  }

  // STEP 2: 3D vapor diffusion + condensation sink
  for (int y = 0; y < yLayers - 1; y++) { 
    for (int x = 1; x < xRes - 1; x++) {
      for (int z = 1; z < zRes - 1; z++) {
        
        float q_c = get_q(x, y, z);
        
        float q_e = (q_c + get_q(x+1, y, z)) * 0.5f;
        float q_w = (q_c + get_q(x-1, y, z)) * 0.5f;
        float flux_x = q_e * (vapor_u[y](x+1, z) - vapor_u[y](x, z)) - 
                       q_w * (vapor_u[y](x, z) - vapor_u[y](x-1, z));

        float q_n = (q_c + get_q(x, y, z+1)) * 0.5f;
        float q_s = (q_c + get_q(x, y, z-1)) * 0.5f;
        float flux_z = q_n * (vapor_u[y](x, z+1) - vapor_u[y](x, z)) - 
                       q_s * (vapor_u[y](x, z) - vapor_u[y](x, z-1));

        float flux_y = 0.0f;
        if (y > 0) {
          float q_u = (q_c + get_q(x, y+1, z)) * 0.5f;
          float q_d = (q_c + get_q(x, y-1, z)) * 0.5f;
          flux_y = q_u * (vapor_u[y+1](x, z) - vapor_u[y](x, z)) - 
                   q_d * (vapor_u[y](x, z) - vapor_u[y-1](x, z));
        } else {
          float q_u = (q_c + get_q(x, y+1, z)) * 0.5f;
          flux_y = q_u * (vapor_u[y+1](x, z) - vapor_u[y](x, z)); 
        }

        // SINK TERM 
        // if at the glass layer, vapor is consumed proportionally to ice growth
        float sink = 0.0f;
        if (y == 0) {
           sink = 0.5f * L_sat * dphi_dt_field(x, z);
        }

        next_u[y](x, z) = vapor_u[y](x, z) + dt * ((flux_x + flux_y + flux_z) / dx2 - sink);
      }
    }
  }

  for (int x = 0; x < xRes; x++)
    for (int z = 0; z < zRes; z++)
      vapor_u[yLayers-1](x,z) = 0.85f; // u_infinity

  // update
  vapor_u = next_u;
  phaseField = next_phase;

  field = phaseField;
}