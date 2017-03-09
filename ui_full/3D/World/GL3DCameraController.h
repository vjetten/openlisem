


#ifndef GL3DCAMERACONTROLLER_H
#define GL3DCAMERACONTROLLER_H

#include <GL3DWidget.h>
#include <GL3DSurface.h>
#include <GL3DCamera.h>

class GL3DSurface;
class GL3DCamera;

class GL3DCameraController
{

public:
    GL3DCameraController(){
    }

    GL3DSurface * m_Surface;

    GL3DCamera * m_Camera;

    bool mode_vertical_up = true;
    bool mode_surface_limit = true;
    bool mode_objects_limit = false;
    bool mode_height_limit = false;
    bool mode_elevation_velocity = true;

    void SetSurface(GL3DSurface * surface);
    void SetMode(int mode);

    void SetCamera(GL3DCamera * c);
    void SetStartPosition();
    void OnUpdate(GL3DWidget * widget,GL3DWorld * world, double dt);
};

#endif // GL3DCAMERACONTROLLER_H
