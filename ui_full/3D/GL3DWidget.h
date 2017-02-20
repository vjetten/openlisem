/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

#include <QtOpenGL/QGLWidget>
#include "QDebug"
#include <QVector>
#include <vector>
#include <QtOpenGL/QtOpenGL>
#include <QOpenGLFunctions_4_0_Core>
#include "csfMap.h"
#include "LisUIoutput.h"
#include "global.h"
#include "3D/World/GL3DWorld.h"
#include "3D/World/GL3DCamera.h"
#include "3D/Graphics/GL3DShaders.h"
#include "3D/Graphics/GL3DTextures.h"
#include "3D/Graphics/GL3DGeometry.h"
#include <QtGlobal>
#include "ui_full/3D/Objects/GL3DObject.h"
#include <QBasicTimer>
#include "3D/Graphics/GL3DMapMath.h"
#include "3D/Graphics/GL3DModels.h"
#include <QMutex>

#ifndef widget3D
#define widget3D

class GL3DMaterial;
class GL3DGeometries;
class GL3DShaders;
class GL3DTextures;
class GL3DModels;
class GL3DCamera;
class GL3DWorld;
class GL3DTexture;
class GL3DGeometry;
class GL3DShader;
class GL3DModel;

#define GL3D_INPUT_NRKEYS 100
#define GL3D_INPUT_NRBUTTONS 5

#define GL3D_MOUSE_LEFT 1
#define GL3D_MOUSE_RIGHT 2
#define GL3D_MOUSE_MIDDLE 2

#define GL3D_KEY_UP 1
#define GL3D_KEY_DOWN 2
#define GL3D_KEY_LEFT 3
#define GL3D_KEY_RIGHT 4
#define GL3D_KEY_W 5
#define GL3D_KEY_A 6
#define GL3D_KEY_S 7
#define GL3D_KEY_D 8
#define GL3D_KEY_CTRL 9
#define GL3D_KEY_SHIFT 10
#define GL3D_KEY_Q 11
#define GL3D_KEY_E 12

#define GL3D_DIR_MODELS "glresources/models/"
#define GL3D_DIR_TEXTURES "glresources/textures/"
#define GL3D_DIR_SHADERS "glresources/shaders/"

#define MV(map,r,c) pcr::isMV(map->data[r][c])
#define INSIDE(map,r, c) (r>=0 && r<map->nrRows() && c>=0 && c<map->nrCols())
#define OUTOFMAP(map,r, c) (r<0 || r >= map->nrRows() || c<0 || c>= map->nrCols())
#define FOR_ROW_COL_MV(mask,map) for(int r = 0; r < map->nrRows(); r++)\
    for (int c = 0; c < map->nrCols(); c++)\
    if(!pcr::isMV(mask->data[r][c]))
#define Drc data[r][c]

#define SAFE_DELETE(x) if(x != 0){delete x;}
class GL3DWidget: public QGLWidget,
        protected QOpenGLFunctions_4_0_Core
{


public:

    bool is_created = false;
    bool ReadyToDraw = false;
    bool m_InitializedRun = false;

    GL3DCamera * m_Camera;
    GL3DWorld * m_World;

    QBasicTimer m_Timer;

    QOpenGLVertexArrayObject m_GLQuadObject;
    QOpenGLVertexArrayObject m_GLQuadObjectChannel;


    QMutex * gl_context_control;
    QOpenGLFunctions_4_0_Core * gl;

    GL3DGeometries * m_Geometries;
    GL3DTextures * m_Textures;
    GL3DShaders * m_Shaders;
    GL3DModels * m_Models;

    QString m_Directory;


    GL3DTexture * m_Texture_Rain;


    GL3DWidget(QWidget *parent = 0):QGLWidget(parent)
    {

    };

    void InitRun(cTMap * copy);
    void SetupShaders();
    void GetData();
    void Update();

    void ClearRun();
    void CloseRun();
    void Close();

    void Onrender();

    void timerEvent(QTimerEvent *e);
    void showEvent(QShowEvent *e);
    void hideEvent(QHideEvent *e);

    void keyPressEvent( QKeyEvent *keyEvent );
    void mouseMoveEvent( QMouseEvent *event );
    void mousePressEvent( QMouseEvent *event );
    void keyReleaseEvent( QKeyEvent *keyEvent );
    void mouseReleaseEvent( QMouseEvent *event );
    void wheelEvent(QWheelEvent *);

    void UseInput();

    //wether or not a key is pressed
    //see GL3D_ enums for key values
    bool KeyPressed[GL3D_INPUT_NRKEYS];
    bool MousePressed[GL3D_INPUT_NRBUTTONS];
    bool KeyPressedT[GL3D_INPUT_NRKEYS];
    bool MousePressedT[GL3D_INPUT_NRBUTTONS];

    //mouse movement since last time that useinput() was called
    bool MouseMovement;        //unit is in screen pixels
    bool MouseMovementTime;    //time in seconds since last time useinput() was called

    int MouseOldX;
    int MouseOldY;
    int MouseX;
    int MouseY;
    int MouseDX;
    int MouseDY;

    int MouseScroll;

    bool KeepMouseMiddle = false;

    int Width = 0;
    int Height = 0;

    qint64 m_Time;
    qint64 m_Time_Start;
    double m_DT;
    double m_Time_s;

    void BindGeometry(QOpenGLVertexArrayObject &object,GL3DShader * s, GL3DGeometry * g);

    void DrawModelObject(GL3DModel * m, GL3DCamera * camera,QMatrix4x4 ModelMatrix);
    void DrawModelGeometryWithMaterial(GL3DGeometry * g,GL3DShader * s,QOpenGLVertexArrayObject * vao,GL3DMaterial * m, GL3DCamera * camera, QMatrix4x4 modelmatrix);
    void DrawModelGeometryWithMaterialMultipleStart(GL3DGeometry * g,GL3DShader * s,QOpenGLVertexArrayObject * vao,GL3DMaterial * m, GL3DCamera * camera);
    void DrawModelGeometryWithMaterialMultiple(GL3DGeometry * g,GL3DShader * s,QOpenGLVertexArrayObject * vao,GL3DMaterial * m, GL3DCamera * camera, QMatrix4x4 modelmatrix);
    void DrawModelGeometryWithMaterialMultipleEnd(GL3DGeometry * g,GL3DShader * s,QOpenGLVertexArrayObject * vao,GL3DMaterial * m, GL3DCamera * camera);

protected:

    void initializeGL();
    void resizeGL(int w, int h);

    void paintGL();

};



class KeyHelper : private QObject {
public:
   static QString keyName(int index) {
      static int keyEnumIndex = staticQtMetaObject.indexOfEnumerator("Key");
      QString name = staticQtMetaObject.enumerator(keyEnumIndex).valueToKey(index);
      if (index >= Qt::Key_Left && index <= Qt::Key_Down) name += " Arrow";
      return name.isEmpty() ? QString() : name.mid(4);
   }
};

#endif
