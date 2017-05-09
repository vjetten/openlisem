#include <GL3DCameraController.h>

void GL3DCameraController::SetSurface(GL3DSurface * surface)
{
    m_Surface = surface;


}

void GL3DCameraController::SetMode(int mode)
{



}

void GL3DCameraController::SetCamera(GL3DCamera * c)
{
    m_Camera = c;

}

void GL3DCameraController::SetStartPosition()
{

   m_Camera->SetPosition(m_Surface->m_Elevation->cellSize() *m_Surface->m_Elevation->nrCols()*0.3,m_Surface->m_ElevationMax + (m_Surface->m_ElevationMax -m_Surface->m_ElevationMin)* 2.0,m_Surface->m_Elevation->cellSize() *m_Surface->m_Elevation->nrRows()*0.3 );
   m_Camera->LookAt(m_Surface->m_Elevation->cellSize() *m_Surface->m_Elevation->nrCols()*0.75,m_Surface->GetElevation(m_Surface->m_Elevation->cellSize() *m_Surface->m_Elevation->nrCols()*0.75,m_Surface->m_Elevation->cellSize() *m_Surface->m_Elevation->nrRows()*0.75 ) ,m_Surface->m_Elevation->cellSize() *m_Surface->m_Elevation->nrRows()*0.75 );

}

void GL3DCameraController::OnUpdate(GL3DWidget * widget,GL3DWorld * world, double dt)
{

    double surf_elevation= m_Surface->GetElevation(m_Camera->m_Position.x(),m_Camera->m_Position.z());
    double elevation = m_Camera->m_Position.y();

    double speed = 15.0 + 2.0 *(elevation - surf_elevation);

    if(widget->KeyPressed[GL3D_KEY_W] == true)
    {
        m_Camera->MoveForward(-dt*speed);
    }
    if(widget->KeyPressed[GL3D_KEY_Q] == true)
    {
        m_Camera->MoveUpward(dt*speed);
    }
    if(widget->KeyPressed[GL3D_KEY_S] == true)
    {
        m_Camera->MoveForward(dt*speed);
    }
    if(widget->KeyPressed[GL3D_KEY_E] == true)
    {
         m_Camera->MoveUpward(-dt*speed);
    }
    if(widget->KeyPressed[GL3D_KEY_A] == true)
    {
        m_Camera->StrafeRight(-dt*speed);
    }
    if(widget->KeyPressed[GL3D_KEY_D] == true)
    {
        m_Camera->StrafeRight(dt*speed);
    }

    if(std::fabs(widget->MouseDX) > 0)
    {
        m_Camera->RotateY(-float(widget->MouseDX));
    }
    if(std::fabs(widget->MouseDY) > 0)
    {
        m_Camera->RotateX(float(widget->MouseDY));
    }

    if(std::fabs(widget->MouseScroll) > 0)
    {
        //m_Camera->Zoom(widget->MouseScroll/25.0);
    }


    surf_elevation= m_Surface->GetElevation(m_Camera->m_Position.x(),m_Camera->m_Position.z());
    elevation = m_Camera->m_Position.y();

    if(elevation < surf_elevation + 0.1 * m_Surface->m_Elevation->cellSize())
    {
        m_Camera->SetPosition(m_Camera->m_Position.x(),surf_elevation + 0.1 * m_Surface->m_Elevation->cellSize(),m_Camera->m_Position.z());

    }

    return;
}
